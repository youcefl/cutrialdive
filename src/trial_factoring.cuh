/*
* Created on 2021.12.23
* Copyright (c) Youcef Lemsafer
*/

#include <cinttypes>
#include <cassert>
#include <algorithm>
#include <string>
#include <iostream>
#include <sstream>
#include <iterator>
#include <vector>
#include <thread>
#include <functional>
#include <omp.h>

#include "device_factors_buffer.hpp"
#include "hgint.hpp"
#include "siever.hpp"
#include "smarandache.hpp"
#include "modular_arithmetic_detail.hpp"
#include "autotests.hpp"


using namespace cutrialdive;


enum class PrecomputeReciprocals {
    no, yes
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
        if(!residue) {
            factorsBuffer.push_factor(numberIndex, divisor);
        }
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


template <uint8_t BatchSize, typename PrimeT, PrecomputeReciprocals precomputeReciprocals>
__global__
void trial_div(uint64_t n0, uint64_t n1,
    uint64_t * firstNumber, uint64_t firstNumberLength,
    data_impl<PrimeT, precomputeReciprocals> primeData,
    typename device_factors_buffer<uint64_t>::device_view_t factorsBuffer
)
{
    static_assert(BatchSize > 0);

    if(n0 >= n1) {
        return;
    }

    int const indexOfFirstPrime = (threadIdx.x + blockIdx.x * blockDim.x) * BatchSize;
    int const stride = blockDim.x * gridDim.x * BatchSize;

    auto primes = primeData.primes;
    auto primesLength = primeData.primes_count;
    for(size_t i = indexOfFirstPrime; i < primesLength; i += stride) {
        uint64_t residue;
        uint64_t divisor;
        if constexpr (BatchSize == 1) {
            divisor = primes[i];
            if constexpr (precomputeReciprocals == PrecomputeReciprocals::yes) {
                residue = modnby1(firstNumber, firstNumberLength, divisor, primeData.reciprocals[i]);
            } else {
                residue = modnby1(firstNumber, firstNumberLength, divisor);
            }
            push_prime_factors_in_batch<BatchSize, PrimeT, precomputeReciprocals>(0, i, primes, primesLength, factorsBuffer, residue, divisor);
        } else {
            divisor = primes[i];
            for(auto j = 1; j < BatchSize; ++j) {
                if (i + j >= primesLength) {
                    break;
                }
                divisor *= primes[i + j];
            }
            residue = modnby1(firstNumber, firstNumberLength, divisor);
            push_prime_factors_in_batch<BatchSize, PrimeT, precomputeReciprocals>(0, i, primes, primesLength, factorsBuffer, residue, divisor);
        }

        /// Propagate residue to next numbers in sequence
        for(auto n = n0 + 1; n < n1; ++n) {
            residue = (__uint128_t(residue) * 100000 + n) % divisor;
            push_prime_factors_in_batch<BatchSize, PrimeT, precomputeReciprocals>(n - n0, i, primes, primesLength, factorsBuffer, residue, divisor);
        }
    }
}

#if 0
// Gets the factors from the device, then repeatedly divides the number by the factors
// to get the cofactor.
// number <- number / (f_1^e_1*f_2^e_2...f_n^e_n)
// @param [in/out] number which was sent to the GPU for factorization
// @return factors (with their exponent)
std::vector<std::pair<uint64_t, uint32_t>> getFactors(HgInt & number, int* deviceFactorsCount, uint64_t * deviceFactors, uint32_t deviceFactorsSize)
{
    std::vector<std::pair<uint64_t, uint32_t>> factors;
    auto const & rawFactors = getFactorsFromDevice(deviceFactorsCount, deviceFactors, deviceFactorsSize);
    for(auto const & f : rawFactors) {
        uint32_t e = 0;
        do {
            mpz_divexact_ui(number.get(), number.get(), f);
            ++e;
        } while (mpz_fdiv_ui(number.get(), f) == 0);
        factors.push_back(std::make_pair(f, e));
    }
    return factors;
}
#endif

auto device_synchronize()
{
    auto cuStatus = cudaDeviceSynchronize();
    if(cuStatus != cudaSuccess) {
        std::cerr << "Error returned by cudaDeviceSynchronize(): " <<
            cudaGetErrorString(cuStatus) << std::endl;
    }
    return cuStatus;
}

void output_usage(std::ostream & out)
{
    out << "Usage:\n    cutrialdive -i <index> --tf-bits <exponent> [--prp-test]" << std::endl;
}


int main(int argc, char** argv)
{
    if(argc == 1) {
        output_usage(std::cout);
        return 0;
    }
    if ((argc == 2) && (std::string(argv[1]) == "-t")) {
        return autotest();
    }
    uint64_t ii = 0;
    int tfBits = 0;
    bool haveToPrpTest{};
    if ((argc >= 3) && (std::string(argv[1]) == "-i")) {
        std::istringstream istri{argv[2]};
        istri >> ii;
        if(argc < 5 || strcmp(argv[3], "--tf-bits")) {
            output_usage(std::cerr);
            return 1;
        }
        std::istringstream istre{argv[4]};
        istre >> tfBits;
        if((tfBits < 2) || (tfBits > 63)) {
            std::cerr << "The value of option --tf-bits must be an integer in range [2, 63]" << std::endl;
            return 1;
        }
        if(argc > 5) {
            if((argc == 6) && !strcmp(argv[5], "--prp-test")) {
                haveToPrpTest = true;
            } else {
                std::cerr << "Unexpected command line argument: `" << argv[5] << "'" << std::endl;
                return 1;
            }
        }
    } else {
        output_usage(std::cerr);
        return 1;
    }
    
    std::cout << "Using up to " << omp_get_max_threads() << " thread(s) on the CPU." << std::endl;

    try {
        auto sm671 = smarandache<HgInt>(ii);
        auto * n = sm671.limbs();
        auto nlen = sm671.size();
        auto sizeInBase10 = sm671.sizeInBase(10);

        std::cout << "Smarandache(" << ii << ") = 0x" << std::hex << n[nlen-1] 
            << "..." << n[0] <<  ", length = " << std::dec << nlen << " 64-bit words, "
            <<  sizeInBase10 << " digits" << std::endl;

        auto sieveUpperBound = uint64_t{1} << tfBits;
        std::cout << "Trial factoring to 2^" << tfBits << std::endl;
        uint64_t * cn;
        auto cuStatus = cudaMalloc(&cn, nlen*sizeof(uint64_t));
        if(cuStatus != cudaSuccess) {
            std::cerr << "Error returned by cudaMalloc(): " << cudaGetErrorString(cuStatus) << std::endl;
        }
        cudaMemcpy(cn, n, nlen*sizeof(uint64_t), cudaMemcpyHostToDevice);


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

        constexpr auto numbersCount = 1;

        cutrialdive::device_factors_buffer<uint64_t> factorsBuffer{ii, numbersCount, 128};

        for (auto k0 = uint64_t{2}, k1 = (std::min)(uint64_t{1} << 25, sieveUpperBound)
                ; k1 <= sieveUpperBound
                ; k0 = k1, k1 += uint64_t{1} << 25
            ) {
            auto sieveStart = std::chrono::high_resolution_clock::now();
            primes.clear();
            cutrialdive::sieve(k0, k1, primes);
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

            if(k1 <= 2642258) {
                trial_div<3, uint64_t><<<256*numSMs, 256>>>(ii, ii + numbersCount, cn, nlen, devicePrimeData.get_data(), factorsBuffer.device_view());
            } else if(k1 <= (1ull << 32)) {
                trial_div<2, uint64_t><<<256*numSMs, 256>>>(ii, ii + numbersCount, cn, nlen, devicePrimeData.get_data(), factorsBuffer.device_view());
            } else {
                trial_div<1, uint64_t><<<256*numSMs, 256>>>(ii, ii + numbersCount, cn, nlen, devicePrimeData.get_data(), factorsBuffer.device_view());
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

        auto foundFactors = factorsBuffer.to_factoring_results<uint64_t, uint32_t>();
        std::cout << foundFactors;

        std::cout << "[Factoring took "
                << std::chrono::duration<double, std::milli>(tfEnd - 
                        tfStart).count() / 1000 << "s ("
                   << "sieve: " << sieve_time.count() << "s)]";
        if(!haveToPrpTest) {
            std::cout << std::endl;
        } else {
            std::cout << " [Primality testing... ";
            std::flush(std::cout);
            auto prpStart = std::chrono::high_resolution_clock::now();
            auto hasPrpCofactor = is_prp(sm671);
            auto prpEnd = std::chrono::high_resolution_clock::now();
            std::cout << std::chrono::duration<double, std::milli>(prpEnd - 
                            prpStart).count() / 1000 << "s]"
                    << std::endl;
            std::cout << "Smarandache(" << ii << ") -> " << (hasPrpCofactor ? "PRP" : "C") << std::endl;
        }
        
    } catch (std::exception const & ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
        return 3;
    }
    return 0;
}

