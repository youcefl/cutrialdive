/*
* Created on 2021.12.23
* Copyright (c) Youcef Lemsafer
*/
#pragma once

#include <cinttypes>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <vector>

#include "device_vector.hpp"
#include "device_factors_buffer.hpp"
#include "prime_data.hpp"
#include "hgint.hpp"
#include "siever.hpp"
#include "modular_arithmetic_detail.hpp"
#include "barrett_mu_types.hpp"
#include "barrett_reciprocals.hpp"
#include "number_sequence.hpp"
#include "number_sequence_helpers.hpp"
#include "progress.hpp"
#include "checkpoint.hpp"
#include "trial_factoring_context.hpp"


namespace cutrialdive {

    template <typename NumberSequenceT, bool HaveInitialValue = InitializeFromValue<NumberSequenceT>>
    struct opt_initial_value
    {
        opt_initial_value(NumberSequenceT, uint64_t) {}
    };

    template <typename NumberSequenceT>
    struct opt_initial_value<NumberSequenceT, true>
    {
        opt_initial_value(NumberSequenceT numSeq, uint64_t n0) {
            auto val = numSeq.value(n0);
            value_at_n0.resize(val.size());
            copy_to_device<uint64_t>(value_at_n0.data(), val.limbs(), val.size());
        }

        device_vector<uint64_t> value_at_n0;
    };

    template <typename NumberSequenceT, typename PrimeT, PrecomputeReciprocals precomputeReciprocals>
    struct device_tf_data : opt_initial_value<NumberSequenceT>
    {
        device_tf_data(NumberSequenceT numSeq,
            uint64_t n0, uint64_t n1,
            data_impl<PrimeT, precomputeReciprocals> primeData,
            typename device_factors_buffer<uint64_t>::device_view_t factorsBuffer
        )
            : opt_initial_value<NumberSequenceT>(numSeq, n0)
            , num_seq(numSeq)
            , n0(n0)
            , n1(n1)
            , prime_data(primeData)
            , factors_buffer(factorsBuffer)
        {}

        NumberSequenceT num_seq;
        uint64_t n0;
        uint64_t n1;
        data_impl<PrimeT, precomputeReciprocals> prime_data;
        typename device_factors_buffer<uint64_t>::device_view_t factors_buffer;

        struct device_view_no_value {
            NumberSequenceT num_seq;
            uint64_t n0;
            uint64_t n1;
            data_impl<PrimeT, precomputeReciprocals> prime_data;
            typename device_factors_buffer<uint64_t>::device_view_t factors_buffer;
        };
        struct device_view_with_value : device_view_no_value {
            uint64_t const* value;
            std::size_t value_length;
        };

        using device_view_t = std::conditional_t<InitializeFromValue<NumberSequenceT>, device_view_with_value, device_view_no_value>;

        device_view_t device_view() const
        {
            if constexpr(InitializeFromValue<NumberSequenceT>) {
                return {num_seq, n0, n1, prime_data, factors_buffer, this->value_at_n0.data(), this->value_at_n0.size()};
            } else {
                return {num_seq, n0, n1, prime_data, factors_buffer};
            }
        }
    };


    __device__ __forceinline__
    void push_prime_factor(
        uint64_t numberIndex,
        typename device_factors_buffer<uint64_t>::device_view_t factorsBuffer,
        uint64_t residue,
        uint64_t divisor
    )
    {
        if(!residue) {
            factorsBuffer.push_factor(numberIndex, divisor);
        }
    }


    template <uint8_t BatchSize, typename PrimeT, PrecomputeReciprocals precomputeReciprocals>
    __device__ __forceinline__
    void push_prime_factors_in_batch(
        uint64_t numberIndex,
        size_t indexOfBatch,
        typename data_impl<PrimeT, precomputeReciprocals>::primes_type primes,
        typename data_impl<PrimeT, precomputeReciprocals>::primes_count_type primesLength,
        typename device_factors_buffer<uint64_t>::device_view_t factorsBuffer,
        uint64_t residue,
        uint64_t divisor
    )
    {
        if constexpr(BatchSize == 1) {
            push_prime_factor(numberIndex, factorsBuffer, residue, divisor);
        } else {
            for (auto j = 0; j < BatchSize; ++j) {
                if (indexOfBatch + j >= primesLength) {
                    break;
                }
                if (! (residue % primes[indexOfBatch + j]) ) {
                    factorsBuffer.push_factor(numberIndex, primes[indexOfBatch + j]);
                }
            }
        }
    }

    template <typename NumberSequenceT, typename PrimeT, PrecomputeReciprocals precomputeReciprocals>
    __global__
    void trial_div_by_2(typename device_tf_data<NumberSequenceT, PrimeT, precomputeReciprocals>::device_view_t data)
    {
        if(data.n0 >= data.n1) {
            return;
        }

        auto residue = value_mod_2(data.num_seq, data.n0);
        push_prime_factor(0, data.factors_buffer, residue, 2);
        /// Propagate residue mod 2 to next numbers in sequence
        for(auto n = data.n0, nEnd = data.n1 - 1; n < nEnd; ++n) {
            residue = next_value_mod_2(data.num_seq, residue, n);
            push_prime_factor(n + 1 - data.n0, data.factors_buffer, residue, 2);
        }
    }

    template <typename NumberSequenceT, uint8_t BatchSize, typename PrimeT, PrecomputeReciprocals precomputeReciprocals>
    __global__
    void trial_div_by_odd_primes(typename device_tf_data<NumberSequenceT, PrimeT, precomputeReciprocals>::device_view_t data)
    {
        static_assert(BatchSize > 0);

        if(data.n0 >= data.n1) {
            return;
        }

        auto primes = data.prime_data.primes;
        auto primesLength = data.prime_data.primes_count;

        int const indexOfFirstPrime = (threadIdx.x + blockIdx.x * blockDim.x) * BatchSize;
        int const stride = blockDim.x * gridDim.x * BatchSize;

        for(size_t i = indexOfFirstPrime; i < primesLength; i += stride) {
            uint64_t divisor = primes[i];
            for(auto j = 1; j < BatchSize; ++j) {
                if (i + j >= primesLength) {
                    break;
                }
                divisor *= primes[i + j];
            }
            uint64_t residue;

            auto barrettMu = compute_barrett_mu_numseq<NumberSequenceT>(divisor);
            if constexpr(InitializeFromValue<NumberSequenceT>) {
                if constexpr ((BatchSize == 1) && (precomputeReciprocals == PrecomputeReciprocals::yes)) {
                    residue = modnby1(data.value, data.value_length, divisor, data.prime_data.reciprocals[i]);
                } else {
                    residue = modnby1(data.value, data.value_length, divisor);
                }
            } else if constexpr(HasValueModMu<NumberSequenceT>) {
                residue = data.num_seq.value_mod_mu(data.n0, divisor, barrettMu);
            } else {
                residue = data.num_seq.value_mod(data.n0, divisor);
            }
            push_prime_factors_in_batch<BatchSize, PrimeT, precomputeReciprocals>(0, i, primes, primesLength, data.factors_buffer, residue, divisor);

            /// Propagate residues to next numbers in sequence
            for(auto n = data.n0, nEnd = data.n1 - 1; n < nEnd; ++n) {
                if constexpr(HasNextValueModMu<NumberSequenceT>) {
                    residue = data.num_seq.next_value_mod_mu(residue, n, divisor, barrettMu);
                } else {
                    residue = data.num_seq.next_value_mod(residue, n, divisor);
                }
                push_prime_factors_in_batch<BatchSize, PrimeT, precomputeReciprocals>(n + 1 - data.n0, i, primes, primesLength, data.factors_buffer, residue, divisor);
            }
        }
    }

    inline
    auto device_synchronize(std::ostream & out)
    {
        auto cuStatus = cudaDeviceSynchronize();
        if(cuStatus != cudaSuccess) {
            out << "Error returned by cudaDeviceSynchronize(): " <<
                cudaGetErrorString(cuStatus) << std::endl;
        }
        return cuStatus;
    }

    inline
    void print_gpu_device_info(int device_id, std::ostream & out)
    {
        cudaDeviceProp props;
        auto cuStatus = cudaGetDeviceProperties(&props, device_id);
        if (cuStatus == cudaSuccess) {
            out << "Using GPU device " << device_id << ": " << props.name
                << " (" << std::setprecision(3) << double(props.totalGlobalMem) / (1024.0 * 1024 * 1024) << " GiB)\n";
        } else {
            out << "Error while getting device info: " << cudaGetErrorString(cuStatus) << std::endl;
        }
    }

    inline
    double delta_now(decltype(std::chrono::high_resolution_clock::now()) start)
    {
        auto now = std::chrono::high_resolution_clock::now();
        if(now <= start) {
            return 0.0;
        }
        return std::chrono::duration<double>(now - start).count();
    }

    template <typename NumberSequenceT>
        requires PureMathSequence<NumberSequenceT>
    inline
    void device_trial_factor(
        trial_factoring_context & ctx,
        NumberSequenceT numSeq
        )
    {
        auto const & opts = ctx.options;
        auto & out = ctx.output_stream;
        auto resumeState = ctx.resume_state;
        auto checkpoint = ctx.checkpoint;
        print_gpu_device_info(0, out);

        int numSMs;
        cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, 0);
        auto tfStart = std::chrono::high_resolution_clock::now();
        double sieveTime{}, computeDataTime{}, copyPrimeDataToHostTime{};
    
        constexpr auto precomputeReciprocals = PrecomputeReciprocals::no;
        prime_data<uint64_t, precomputeReciprocals> primeData;
        primeData.reserve(1ull << 22);
        device_prime_data<uint64_t, precomputeReciprocals> devicePrimeData;

        device_factors_buffer<uint64_t> factorsBuffer = resumeState 
                ? device_factors_buffer<uint64_t>{resumeState->factors_buffer.get()}
                : device_factors_buffer<uint64_t>{opts.n0, opts.n1 - opts.n0, opts.max_factors_per_number};

        if(resumeState) {
            out << "Resuming at the smallest prime > " << resumeState->last_processed_prime << "." << std::endl;
        }

        device_tf_data<NumberSequenceT, uint64_t, precomputeReciprocals> tf_data{
            numSeq, opts.n0, opts.n1, devicePrimeData.get_data(), factorsBuffer.device_view()};

        auto progressHandler = opts.is_progress_enabled ? std::make_unique<progress>(opts.f1, ctx.runtime_options.progress_period, out) 
                                                        : std::unique_ptr<progress>{};
        auto updateProgress = progressHandler
            ? std::function<void(uint64_t, uint64_t)>{[&progressHandler](auto segUpperBound, auto lastPrimeInSeg) {
                progressHandler->update(segUpperBound, lastPrimeInSeg); 
              }}
            : std::function<void(uint64_t, uint64_t)>{[](auto, auto) {}};

        // Special treatment for 2 if needed
        auto f0 = resumeState ? resumeState->last_processed_prime + 1 : opts.f0;
        if(f0 <= 2 && opts.f1 > 2) {
            f0 = 3;
            trial_div_by_2<NumberSequenceT, uint64_t, precomputeReciprocals><<<1, 1>>>(tf_data.device_view());
        }

        constexpr auto segmentLen = uint64_t{1} << 26;
        auto const f1 = opts.f1;
        auto lastPrimeOfPreviousSeg = uint64_t{};
        siever primeGen{ctx.runtime_options.threads_count};

        out << "Sieving for primes using ";
        if(primeGen.threads_count() == 1) {
            out << "one thread." << std::endl;
        } else {
            out << "up to " << primeGen.threads_count() << " threads." << std::endl;
        }
        std::vector<std::span<uint64_t>> primes;
        primes.reserve(primeGen.threads_count());

        // Main loop for odd primes only
        for (auto fx = f0, fy = (std::min)(f0 + segmentLen, f1);
            fx < f1;
            fx = fy, fy = (std::min)(fy + segmentLen, f1)
        ) {
            auto sieveStart = std::chrono::high_resolution_clock::now();
            primeGen.sieve(fx, fy, primes);
            sieveTime += delta_now(sieveStart);

            if constexpr(precomputeReciprocals == PrecomputeReciprocals::yes) {
                auto computeDataStart = std::chrono::high_resolution_clock::now();
                compute_prime_data(primes, primeData);
                computeDataTime += delta_now(computeDataStart);
            }

            device_synchronize(out);

            if(lastPrimeOfPreviousSeg) {
                updateProgress(fx, lastPrimeOfPreviousSeg);
                if(checkpoint && checkpoint->due()) {
                    checkpoint->write(engine_state{
                        lastPrimeOfPreviousSeg,
                        factors_buffer_holder{factorsBuffer.to_host_factors_buffer<uint64_t>()}
                    });
                }
            }

            auto copyPrimeDataToHostStart = std::chrono::high_resolution_clock::now();
            devicePrimeData.copy_from_host(primes, primeData);
            copyPrimeDataToHostTime += delta_now(copyPrimeDataToHostStart);
            tf_data.prime_data = devicePrimeData.get_data();

            if(fy <= 2642258) {
                trial_div_by_odd_primes<NumberSequenceT, 3, uint64_t, precomputeReciprocals><<<256*numSMs, 256>>>(tf_data.device_view());
            } else if(fy <= (1ull << 32)) {
                trial_div_by_odd_primes<NumberSequenceT, 2, uint64_t, precomputeReciprocals><<<256*numSMs, 256>>>(tf_data.device_view());
            } else {
                trial_div_by_odd_primes<NumberSequenceT, 1, uint64_t, precomputeReciprocals><<<256*numSMs, 256>>>(tf_data.device_view());
            }
            auto err = cudaGetLastError();
            if(err != cudaSuccess) {
                out << "Error executing trial division kernel" << std::endl;
                out << "Launch error: " << cudaGetErrorString(err) << std::endl;
            }
            lastPrimeOfPreviousSeg = !primes.empty() ? primes.back().back() : 0;
        }
        device_synchronize(out);
        if(progressHandler) {
            progressHandler->end();
        }
        if(checkpoint) {
            checkpoint->end();
        }

        auto tfEnd = std::chrono::high_resolution_clock::now();
        

        if constexpr(precomputeReciprocals == PrecomputeReciprocals::yes) {
            out << "Computing prime data on host took " << computeDataTime << "s (cumulated time)" << std::endl;
        }
        out << "Copying prime data to device took " << copyPrimeDataToHostTime << "s (cumulated time)" << std::endl;

        ctx.results = factorsBuffer.to_factoring_results<uint64_t, uint32_t>();

        out << "[Factoring took "
                << std::chrono::duration<double, std::milli>(tfEnd - 
                        tfStart).count() / 1000 << "s ("
                   << "sieve: " << sieveTime << "s)]" << std::endl;

    }

} // namespace cutrialdive

