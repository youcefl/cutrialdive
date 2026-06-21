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
#include "logger.hpp"


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

    struct kernel_launch_config
    {
        uint8_t batch_size;
        int grid_size;
        int block_size;
        double max_occupancy;
    };

    template <typename NumberSequenceT, typename PrimeT, PrecomputeReciprocals precomputeReciprocals>
        requires PureMathSequence<NumberSequenceT>
    inline
    std::vector<kernel_launch_config> get_occupancy_info(cudaDeviceProp & prop, tf_runtime_options const & rtOpts)
    {
        if(rtOpts.grid_size.has_value() && rtOpts.block_size.has_value()) {
            std::vector<kernel_launch_config> result;
            for(int batchSize = 1; batchSize <= 3; ++batchSize) {
                result.emplace_back(batchSize, *rtOpts.grid_size, *rtOpts.block_size);
            }
            return result;
        }
        auto kernels = {
            trial_div_by_odd_primes<NumberSequenceT, 1, uint64_t, precomputeReciprocals>,
            trial_div_by_odd_primes<NumberSequenceT, 2, uint64_t, precomputeReciprocals>,
            trial_div_by_odd_primes<NumberSequenceT, 3, uint64_t, precomputeReciprocals>
        };
        std::vector<kernel_launch_config> result;
        result.reserve(3);

        uint8_t batchSize{1};
        for(auto & ker : kernels) {
            int blockSize = rtOpts.block_size.value_or([&](){
                int minGridSize{};
                int tmpBlockSize{};
                auto err = cudaOccupancyMaxPotentialBlockSize(&minGridSize, &tmpBlockSize, ker);
                if(err != cudaSuccess) {
                    tmpBlockSize = 256;
                    std::cout << "Error while getting block size:" << std::endl;
                    std::cout << cudaGetErrorString(err) << std::endl;
                    std::cout << "Using default value " << tmpBlockSize << std::endl;
                }
                return tmpBlockSize;
            }());
            
            int gridSize = rtOpts.grid_size.value_or([&](){
                int maxBlocksPerSM;
                auto err = cudaOccupancyMaxActiveBlocksPerMultiprocessor(
                            &maxBlocksPerSM, ker, blockSize, 0);
                if(err != cudaSuccess) {
                    maxBlocksPerSM = 64;
                    std::cout << "Error while getting maximum blocks per SM:" << std::endl;
                    std::cout << cudaGetErrorString(err) << std::endl;
                    std::cout << "Using default value " << maxBlocksPerSM << std::endl;
                }
                return prop.multiProcessorCount * maxBlocksPerSM * 4;             
            }());
    
            result.emplace_back(batchSize, gridSize, blockSize);

            CTDLOG_DEBUG() << "Kernel launch config for batch size " << uint32_t(batchSize) 
                << ": <<<" << gridSize << ", " << blockSize << ">>>" << std::endl;

            ++batchSize;
        }
        return result;
    }


    template <typename NumberSequenceT>
        requires PureMathSequence<NumberSequenceT>
    inline
    void device_trial_factor(
        trial_factoring_context & ctx,
        NumberSequenceT numSeq
        )
    {
        int const deviceId = ctx.runtime_options.device_id;
        CTDLOG_DEBUG() << "About to perform trial factoring on device " << deviceId << "...\n" << ctx.options << std::endl;
        auto const & opts = ctx.options;
        auto & out = ctx.output_stream;
        auto resumeState = ctx.resume_state;
        auto checkpoint = ctx.checkpoint;

        CTDLOG_DIAG() << "Calling cudaSetDevice(" << deviceId << ")" << std::endl;
        auto cudaErr = cudaSetDevice(deviceId);
        if(cudaErr != cudaSuccess) {
            out << "Error selecting device " << deviceId << ": "
                << cudaGetErrorString(cudaErr) << std::endl;
            return;
        }
        print_gpu_device_info(deviceId, out);
        cudaDeviceProp prop;
        cudaErr = cudaGetDeviceProperties(&prop, deviceId);
        if(cudaErr != cudaSuccess) {
            out << "Error getting device properties: "
                << cudaGetErrorString(cudaErr) << std::endl;
            return;
        }
        CTDLOG_VERBOSE() << "Number of SMs: " << prop.multiProcessorCount << std::endl;
        if(ctx.runtime_options.grid_size.has_value()) {
            CTDLOG_INFO() << "Using non-default value " << *ctx.runtime_options.grid_size << " for grid size." << std::endl;
        }
        if(ctx.runtime_options.block_size.has_value()) {
            CTDLOG_INFO() << "Using non-default value " << *ctx.runtime_options.block_size << " for block size." << std::endl;
        }


        auto tfStart = std::chrono::high_resolution_clock::now();
        double sieveTime{}, computeDataTime{}, copyPrimeDataToHostTime{};
    
        constexpr auto precomputeReciprocals = PrecomputeReciprocals::no;
        auto kernelConfigs = get_occupancy_info<NumberSequenceT, uint64_t, precomputeReciprocals>(
            prop,
            ctx.runtime_options
          );
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

        auto f0 = resumeState ? resumeState->last_processed_prime + 1 : opts.f0;
        auto progressHandler = ctx.runtime_options.is_progress_enabled
                ? std::make_unique<progress>(f0, opts.f1, ctx.runtime_options.progress_period, out) 
                : std::unique_ptr<progress>{};
        auto updateProgress = progressHandler
            ? std::function<void(uint64_t, uint64_t)>{[&progressHandler](auto segUpperBound, auto lastPrimeInSeg) {
                progressHandler->update(segUpperBound, lastPrimeInSeg); 
              }}
            : std::function<void(uint64_t, uint64_t)>{[](auto, auto) {}};

        // Special treatment for 2 if needed
        if(f0 <= 2 && opts.f1 > 2) {
            f0 = 3;
            trial_div_by_2<NumberSequenceT, uint64_t, precomputeReciprocals><<<1, 1>>>(tf_data.device_view());
        }

        auto const segmentLen = ctx.runtime_options.segment_length;
        CTDLOG_VERBOSE() << "Length of sieve segment: " << segmentLen << std::endl;
        auto const f1 = opts.f1;
        auto lastPrimeOfPreviousSeg = uint64_t{};
        siever primeGen{ctx.runtime_options.common_options.threads_count};

        out << "Sieving for primes using ";
        if(primeGen.threads_count() == 1) {
            out << "one thread." << std::endl;
        } else {
            out << "up to " << primeGen.threads_count() << " threads." << std::endl;
        }
        std::vector<std::span<uint64_t>> primes;
        primes.reserve(primeGen.threads_count());

        // Main loop for odd primes only
        decltype(opts.f1) lastSegEnd{};
        for (auto fx = f0, fy = (std::min)(f0 + segmentLen, f1);
            fx < f1;
            fx = fy, fy = (std::min)(fy + segmentLen, f1)
        ) {
            CTDLOG_DIAG() << "Sieving segment [" << fx << ", " << fy << "[" << "\n";
            lastSegEnd = fy;
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
                trial_div_by_odd_primes<NumberSequenceT, 3, uint64_t, precomputeReciprocals>
                    <<<kernelConfigs[2].grid_size, kernelConfigs[2].block_size>>>(tf_data.device_view());
            } else if(fy <= (1ull << 32)) {
                trial_div_by_odd_primes<NumberSequenceT, 2, uint64_t, precomputeReciprocals>
                    <<<kernelConfigs[1].grid_size, kernelConfigs[1].block_size>>>(tf_data.device_view());
            } else {
                trial_div_by_odd_primes<NumberSequenceT, 1, uint64_t, precomputeReciprocals>
                    <<<kernelConfigs[0].grid_size, kernelConfigs[0].block_size>>>(tf_data.device_view());
            }
            auto err = cudaGetLastError();
            if(err != cudaSuccess) {
                out << "Error executing trial division kernel" << std::endl;
                out << "Launch error: " << cudaGetErrorString(err) << std::endl;
            }
            lastPrimeOfPreviousSeg = !primes.empty() ? primes.back().back() : 0;
        }
        device_synchronize(out);
        if(lastSegEnd && lastPrimeOfPreviousSeg) {
            updateProgress(lastSegEnd, lastPrimeOfPreviousSeg);
        }
        if(progressHandler) {
            progressHandler->end();
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

