/*
* Created on 2021.12.23
* Copyright (c) Youcef Lemsafer
*/
#pragma once
#include <cinttypes>
#include <iostream>
#include <vector>

#include "device_vector.hpp"
#include "device_factors_buffer.hpp"
#include "hgint.hpp"
#include "siever.hpp"
#include "modular_arithmetic_detail.hpp"
#include "barrett_mu_types.hpp"
#include "barrett_reciprocals.hpp"
#include "number_sequence.hpp"


namespace cutrialdive {

    enum class PrecomputeReciprocals
    {
        no,
        yes
    };

    template <typename PrimeT, PrecomputeReciprocals havePrecomputedReciprocals>
    class prime_data;

    template <typename PrimeT>
    class prime_data<PrimeT, PrecomputeReciprocals::no>
    {
    public:
        uint64_t size() const;
        void reserve(size_t capacity);
    };
    template <typename PrimeT>
    inline
    uint64_t prime_data<PrimeT, PrecomputeReciprocals::no>::size() const
    {
        return 0;
    }
    template <typename PrimeT>
    inline
    void prime_data<PrimeT, PrecomputeReciprocals::no>::reserve(std::size_t )
    {
    }


    template <typename PrimeT>
    class prime_data<PrimeT, PrecomputeReciprocals::yes>
    {
    public:
        uint64_t size() const;
        void reserve(size_t capacity);
        std::vector<PrimeT> const & reciprocals() const;
        std::vector<PrimeT> & reciprocals();
    private:
        std::vector<PrimeT> reciprocals_;
    };
    template <typename PrimeT>
    inline
    uint64_t prime_data<PrimeT, PrecomputeReciprocals::yes>::size() const
    {
        return reciprocals_.size();
    }
    template <typename PrimeT>
    inline
    void prime_data<PrimeT, PrecomputeReciprocals::yes>::reserve(std::size_t capacity)
    {
        reciprocals_.reserve(capacity);
    }
    template <typename PrimeT>
    inline
    std::vector<PrimeT> const & prime_data<PrimeT, PrecomputeReciprocals::yes>::reciprocals() const
    {
        return reciprocals_;
    }
    template <typename PrimeT>
    inline
    std::vector<PrimeT> & prime_data<PrimeT, PrecomputeReciprocals::yes>::reciprocals()
    {
        return reciprocals_;
    }

    template <typename PrimeT, PrecomputeReciprocals precomputeReciprocals>
    inline void compute_prime_data(std::vector<PrimeT> const & primes, prime_data<PrimeT, precomputeReciprocals> & data)
    {
        if constexpr (precomputeReciprocals == PrecomputeReciprocals::yes) {
            cutrialdive::compute_reciprocals(primes, data.reciprocals());
        }
    }

    template <typename PrimeT, PrecomputeReciprocals precomputeReciprocals>
    struct data_impl;
    template <typename PrimeT>
    struct data_impl<PrimeT, PrecomputeReciprocals::yes> {
        using primes_type = PrimeT*;
        using primes_count_type = uint64_t;

        PrimeT * primes;
        PrimeT * reciprocals;
        uint64_t primes_count;
    };
    template <typename PrimeT>
    struct data_impl<PrimeT, PrecomputeReciprocals::no> {
        using primes_type = PrimeT*;
        using primes_count_type = uint64_t;

        primes_type primes;
        primes_count_type primes_count;
    };

    template <typename PrimeT, PrecomputeReciprocals precomputeReciprocals>
    class device_prime_data {
    public:
        device_prime_data() = default;
        device_prime_data(device_prime_data&&) = default;
        device_prime_data(device_prime_data const &) = delete;
        device_prime_data& operator=(device_prime_data const &) = delete;

        data_impl<PrimeT, precomputeReciprocals> const & get_data() const;
        uint64_t size() const;
        void resize(uint64_t newSize);
        void copy_from_host(
                std::vector<PrimeT> const & primes,
                prime_data<PrimeT, precomputeReciprocals> const & data
            );
    private:
        data_impl<PrimeT, precomputeReciprocals> data_{};
    };




    template <typename PrimeT, PrecomputeReciprocals precomputeReciprocals>
    inline
    data_impl<PrimeT, precomputeReciprocals> const & device_prime_data<PrimeT, precomputeReciprocals>::get_data() const
    {
        return data_;   
    }


    template <typename PrimeT, PrecomputeReciprocals precomputeReciprocals>
    inline uint64_t device_prime_data<PrimeT, precomputeReciprocals>::size() const
    {
    //    std::cout << "Size of device_prime_data instance: " << data_.primes_count << std::endl;
        return data_.primes_count;
    }

    template <typename PrimeT, PrecomputeReciprocals precomputeReciprocals>
    inline void device_prime_data<PrimeT, precomputeReciprocals>::resize(uint64_t newSize)
    {
    //    std::cout << "Resizing data in device_prime_data instance: {" 
    //        << data_.primes << ", " << data_.normalized
    //        << ", " << data_.reciprocals << "}" << std::endl;
        if(newSize <= data_.primes_count) {
            data_.primes_count = newSize;
            return;
        }
        data_.primes_count = newSize;
        free_device_ptr(data_.primes);
        data_.primes = alloc_device_ptr<PrimeT>(data_.primes_count);
        if constexpr (precomputeReciprocals == PrecomputeReciprocals::yes) {
            free_device_ptr(data_.reciprocals);
            data_.reciprocals = alloc_device_ptr<PrimeT>(data_.primes_count);
        }
    }

    template <typename PrimeT, PrecomputeReciprocals precomputeReciprocals>
    inline void device_prime_data<PrimeT, precomputeReciprocals>::copy_from_host(
        std::vector<PrimeT> const & primes,
        prime_data<PrimeT, precomputeReciprocals> const & data
        )
    {
        if constexpr (precomputeReciprocals == PrecomputeReciprocals::yes) {
            // Must have same number of primes and additional data
            // (because there must be a one to one correspondance between
            // i-th prime and associated data)
            assert(primes.size() == data.size());
        }
        auto const primes_size = primes.size();
        resize(primes_size);
        if(!size()) {
            return;
        }
        cudaMemcpy(data_.primes, &primes[0], sizeof(*&primes[0]) * primes_size, cudaMemcpyHostToDevice);
        if constexpr (precomputeReciprocals == PrecomputeReciprocals::yes) {
            cudaMemcpy(data_.reciprocals, &data.reciprocals()[0],
                    sizeof(*&data.reciprocals()[0]) * data.reciprocals().size(),
                    cudaMemcpyHostToDevice
                );
        }
    }


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

    template <typename NumberSequenceT, typename PrimeT>
    __global__
    void trial_div_by_2(uint64_t n0, uint64_t n1, typename device_factors_buffer<uint64_t>::device_view_t factorsBuffer)
    {
        if(n0 >= n1) {
            return;
        }

        auto residue = NumberSequenceT::value_mod_2(n0);
        push_prime_factor(0, factorsBuffer, residue, 2);
        /// Propagate residue mod 2 to next numbers in sequence
        for(auto n = n0, nEnd = n1 - 1; n < nEnd; ++n) {
            residue = NumberSequenceT::next_value_mod_2(residue, n);
            push_prime_factor(n + 1 - n0, factorsBuffer, residue, 2);
        }
    }

    template <typename NumberSequenceT, uint8_t BatchSize, typename PrimeT, PrecomputeReciprocals precomputeReciprocals>
    __global__
    void trial_div_by_odd_primes(uint64_t n0, uint64_t n1,
        uint64_t * firstNumber, uint64_t firstNumberLength,
        data_impl<PrimeT, precomputeReciprocals> primeData,
        typename device_factors_buffer<uint64_t>::device_view_t factorsBuffer
    )
    {
        static_assert(BatchSize > 0);

        if(n0 >= n1) {
            return;
        }

        auto primes = primeData.primes;
        auto primesLength = primeData.primes_count;

        int const indexOfFirstPrime = (threadIdx.x + blockIdx.x * blockDim.x) * BatchSize;
        int const stride = blockDim.x * gridDim.x * BatchSize;


        for(size_t i = indexOfFirstPrime; i < primesLength; i += stride) {
            uint64_t residue;
            uint64_t divisor = primes[i];
            if constexpr (BatchSize == 1) {
                if constexpr (precomputeReciprocals == PrecomputeReciprocals::yes) {
                    residue = modnby1(firstNumber, firstNumberLength, divisor, primeData.reciprocals[i]);
                } else {
                    residue = modnby1(firstNumber, firstNumberLength, divisor);
                }
                push_prime_factor(0, factorsBuffer, residue, divisor);
            } else {
                for(auto j = 1; j < BatchSize; ++j) {
                    if (i + j >= primesLength) {
                        break;
                    }
                    divisor *= primes[i + j];
                }
                residue = modnby1(firstNumber, firstNumberLength, divisor);
                push_prime_factors_in_batch<BatchSize, PrimeT, precomputeReciprocals>(0, i, primes, primesLength, factorsBuffer, residue, divisor);
            }

            /// Propagate residues to next numbers in sequence
            auto barrettMu = compute_barrett_mu<NumberSequenceT>(divisor);
            for(auto n = n0, nEnd = n1 - 1; n < nEnd; ++n) {
                if constexpr(std::is_same_v<decltype(barrettMu), no_barrett_t>) {
                    residue = NumberSequenceT::next_value_mod(residue, n, divisor);
                } else {
                    residue = NumberSequenceT::next_value_mod_mu(residue, n, divisor, barrettMu);
                }
                push_prime_factors_in_batch<BatchSize, PrimeT, precomputeReciprocals>(n + 1 - n0, i, primes, primesLength, factorsBuffer, residue, divisor);
            }
        }
    }

    inline
    auto device_synchronize()
    {
        auto cuStatus = cudaDeviceSynchronize();
        if(cuStatus != cudaSuccess) {
            std::cerr << "Error returned by cudaDeviceSynchronize(): " <<
                cudaGetErrorString(cuStatus) << std::endl;
        }
        return cuStatus;
    }

    inline
    void print_gpu_device_info(int device_id = 0) {
        cudaDeviceProp props;
        auto cuStatus = cudaGetDeviceProperties(&props, device_id);
        if (cuStatus == cudaSuccess) {
            printf("GPU device %d: %s (%.1f GiB)\n", 
                device_id, 
                props.name, 
                double(props.totalGlobalMem) / (1024.0 * 1024 * 1024));
        } else {
            std::cerr << "Error while getting device info: " << cudaGetErrorString(cuStatus) << std::endl;
        }
    }

    template <typename NumberSequenceT>
    inline
    void device_trial_factor(
        trial_factoring_options const & opts,
        factoring_results<uint64_t, uint32_t> & results
        )
    {

        print_gpu_device_info();
        auto sn0 = NumberSequenceT::value(opts.n0);
        auto * sn0Limbs = sn0.limbs();
        auto sn0Len = sn0.size();

        uint64_t * cuN;
        auto const cnBytesCount = sn0Len * sizeof(*sn0Limbs);
        auto cuStatus = cudaMalloc(&cuN, cnBytesCount);
        if(cuStatus != cudaSuccess) {
            std::cerr << "Error returned by cudaMalloc(): " << cudaGetErrorString(cuStatus) << std::endl;
        }
        cudaMemcpy(cuN, sn0Limbs, cnBytesCount, cudaMemcpyHostToDevice);

        int numSMs;
        cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, 0);
        auto tfStart = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> sieve_time, computeDataTime, copyPrimeDataToHostTime;
    
        constexpr auto precomputeReciprocals = PrecomputeReciprocals::no;
        std::vector<uint64_t> primes;
        primes.reserve(1ull << 22);
        prime_data<uint64_t, precomputeReciprocals> primeData;
        primeData.reserve(1ull << 22);
        device_prime_data<uint64_t, precomputeReciprocals> devicePrimeData;

        cutrialdive::device_factors_buffer<uint64_t> factorsBuffer{opts.n0, opts.n1 - opts.n0, opts.max_factors_per_number};

        // Special treatment for 2 if needed
        auto f0 = opts.f0;
        if(f0 <= 2 && opts.f1 > 2) {
            f0 = 3;
            trial_div_by_2<NumberSequenceT, uint64_t><<<1, 1>>>(opts.n0, opts.n1, factorsBuffer.device_view());
        }

        constexpr auto segmentLen = uint64_t{1} << 25;
        auto const f1 = opts.f1;

        // Main loop for odd primes only
        for (auto fx = f0, fy = (std::min)(f0 + segmentLen, f1);
            fx < f1;
            fx = fy, fy = (std::min)(fy + segmentLen, f1)
        ) {
            auto sieveStart = std::chrono::high_resolution_clock::now();
            primes.clear();
            sieve(fx, fy, primes);
            sieve_time += std::chrono::high_resolution_clock::now() - sieveStart;

            if constexpr(precomputeReciprocals == PrecomputeReciprocals::yes) {
                auto computeDataStart = std::chrono::high_resolution_clock::now();
                compute_prime_data(primes, primeData);
                computeDataTime += std::chrono::high_resolution_clock::now() - computeDataStart;
            }

            device_synchronize();

            auto copyPrimeDataToHostStart = std::chrono::high_resolution_clock::now();
            devicePrimeData.copy_from_host(primes, primeData);
            copyPrimeDataToHostTime += std::chrono::high_resolution_clock::now() - copyPrimeDataToHostStart;

            if(fy <= 2642258) {
                trial_div_by_odd_primes<NumberSequenceT, 3, uint64_t><<<256*numSMs, 256>>>(opts.n0, opts.n1, cuN, sn0Len, devicePrimeData.get_data(), factorsBuffer.device_view());
            } else if(fy <= (1ull << 32)) {
                trial_div_by_odd_primes<NumberSequenceT, 2, uint64_t><<<256*numSMs, 256>>>(opts.n0, opts.n1, cuN, sn0Len, devicePrimeData.get_data(), factorsBuffer.device_view());
            } else {
                trial_div_by_odd_primes<NumberSequenceT, 1, uint64_t><<<256*numSMs, 256>>>(opts.n0, opts.n1, cuN, sn0Len, devicePrimeData.get_data(), factorsBuffer.device_view());
            }
            auto err = cudaGetLastError();
            if(err != cudaSuccess) {
                std::cout << "Error executing kernel" << std::endl;
                printf("Launch error: %s\n", cudaGetErrorString(err));
            }
        }
        device_synchronize();

        auto tfEnd = std::chrono::high_resolution_clock::now();
        

        if constexpr(precomputeReciprocals == PrecomputeReciprocals::yes) {
            std::cout << "Computing prime data on host took " << computeDataTime.count() << "s (cumulated time)" << std::endl;
        }
        std::cout << "Copying prime data to device took " << copyPrimeDataToHostTime.count() << "s (cumulated time)" << std::endl;

        results = factorsBuffer.to_factoring_results<uint64_t, uint32_t>();

        std::cout << "[Factoring took "
                << std::chrono::duration<double, std::milli>(tfEnd - 
                        tfStart).count() / 1000 << "s ("
                   << "sieve: " << sieve_time.count() << "s)]" << std::endl;

    }

} // namespace cutrialdive

