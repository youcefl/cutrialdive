/*
* Youcef Lemsafer
* 2021.12.23
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

#include "hgint.hpp"
#include "siever.hpp"
#include "smarandache.hpp"
#include "modular_arithmetic_detail.hpp"
#include "autotests.hpp"


using namespace cutrialdive;

template <typename PtrT>
void free_device_ptr(PtrT *& ptr)
{
    if(!ptr) {
        return;
    }
    cudaFree(ptr);
    ptr = nullptr;
}

template <typename PtrT>
PtrT* alloc_device_ptr(uint64_t size)
{
    PtrT * ptr = nullptr;
    cudaMalloc(&ptr, sizeof(*ptr) * size);
    return ptr;
}


template <typename PrimeT>
class prime_data
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
uint64_t prime_data<PrimeT>::size() const
{
    return reciprocals_.size();
}
template <typename PrimeT>
inline
void prime_data<PrimeT>::reserve(std::size_t capacity)
{
    reciprocals_.reserve(capacity);
}
template <typename PrimeT>
inline
std::vector<PrimeT> const & prime_data<PrimeT>::reciprocals() const
{
    return reciprocals_;
}
template <typename PrimeT>
inline
std::vector<PrimeT> & prime_data<PrimeT>::reciprocals()
{
    return reciprocals_;
}

template <typename PrimeT>
inline void compute_prime_data(std::vector<PrimeT> const & primes, prime_data<PrimeT> & data)
{
    cutrialdive::compute_reciprocals(primes, data.reciprocals());
}

template <typename PrimeT>
class device_prime_data {
public:
    struct data {
        PrimeT * primes;
        PrimeT * reciprocals;
        uint64_t primes_count;
    };
    device_prime_data() = default;
    device_prime_data(device_prime_data&&) = default;
    device_prime_data(device_prime_data const &) = delete;
    device_prime_data& operator=(device_prime_data const &) = delete;

    data const & get_data() const;
    uint64_t size() const;
    void resize(uint64_t newSize);
    void copy_from_host(
            std::vector<PrimeT> const & primes,
            prime_data<PrimeT> const & data
        );
private:
    data data_{};
};

template <typename PrimeT>
inline
typename device_prime_data<PrimeT>::data const & device_prime_data<PrimeT>::get_data() const
{
    return data_;   
}


template <typename PrimeT>
inline uint64_t device_prime_data<PrimeT>::size() const
{
//    std::cout << "Size of device_prime_data instance: " << data_.primes_count << std::endl;
    return data_.primes_count;
}

template <typename PrimeT>
inline void device_prime_data<PrimeT>::resize(uint64_t newSize)
{
//    std::cout << "Resizing data in device_prime_data instance: {" 
//        << data_.primes << ", " << data_.normalized
//        << ", " << data_.reciprocals << "}" << std::endl;
    free_device_ptr(data_.primes);
    free_device_ptr(data_.reciprocals);
    data_.primes_count = newSize;
    data_.primes = alloc_device_ptr<PrimeT>(data_.primes_count);
    data_.reciprocals = alloc_device_ptr<PrimeT>(data_.primes_count);
}

template <typename PrimeT>
inline void device_prime_data<PrimeT>::copy_from_host(
    std::vector<PrimeT> const & primes,
    prime_data<PrimeT> const & data
    )
{
    assert(primes.size() == data.size());
    auto const primes_size = primes.size();
    if(primes_size > size()) {
        resize(primes_size);
    }
    if(!size()) {
        return;
    }
    cudaMemcpy(data_.primes, &primes[0], sizeof(*&primes[0]) * primes_size, cudaMemcpyHostToDevice);
    cudaMemcpy(data_.reciprocals, &data.reciprocals()[0], sizeof(*&data.reciprocals()[0]) * primes_size, cudaMemcpyHostToDevice);
}

__device__
void pushFactor(uint64_t ff, int* deviceFactorsCount, uint64_t * deviceFactors, uint32_t deviceFactorsSize)
{
//    __syncthreads();
    auto factorIdx = atomicAdd(deviceFactorsCount, 1);
    deviceFactors[factorIdx] = ff;
//    printf("Device factors count = %d\n", *deviceFactorsCount);
//    __syncthreads();
}

template <uint8_t BatchSize, typename PrimeT>
__global__
void trial_div(uint64_t * n, uint64_t nlen, typename device_prime_data<PrimeT>::data primeData, int* deviceFactorsCount, uint64_t * deviceFactors, uint32_t deviceFactorsSize)
{
    static_assert(BatchSize > 0);
    int index = (threadIdx.x + blockIdx.x * blockDim.x) * BatchSize;
    int stride = blockDim.x * gridDim.x * BatchSize;

    auto p = primeData.primes;
    auto plen = primeData.primes_count;
    for(size_t i = index; i < plen; i += stride) {
        if constexpr (BatchSize == 1) {
            if(!modnby1(n, nlen, p[i], primeData.reciprocals[i])) {
                pushFactor(p[i], deviceFactorsCount, deviceFactors, deviceFactorsSize);
            }
        } else {
            uint64_t pp = p[i];
            for(auto j = 1; j < BatchSize; ++j) {
                if (i + j >= plen) {
                    break;
                }
                pp *= p[i + j];
            }
            auto mod = modnby1(n, nlen, pp);
            for (auto j = 0; j < BatchSize; ++j) {
                if (i + j >= plen) {
                    break;
                }
                if (! (mod % p[i + j]) ) {
    //    printf("Adding factor %llu\n", p[i+j]);
                    pushFactor(p[i + j], deviceFactorsCount, deviceFactors, deviceFactorsSize);
                }
            }
        }
    }
}

std::vector<uint64_t> getFactorsFromDevice(int* deviceFactorsCount, uint64_t * deviceFactors, uint32_t deviceFactorsSize)
{
    int factorsCount = 0;
    cudaMemcpy(&factorsCount, deviceFactorsCount, sizeof(int), cudaMemcpyDeviceToHost);
//    std::cout << "factorsCount = " << factorsCount << std::endl;
    auto ret = std::vector<uint64_t>(std::size_t(factorsCount), uint64_t(0));
//    std::cout << "Size of ret = " << ret.size() << std::endl;
    cudaMemcpy(&ret[0], deviceFactors, sizeof(uint64_t)*factorsCount, cudaMemcpyDeviceToHost);
    std::sort(std::begin(ret), std::end(ret));
    return ret;
}







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

void alloc_device_factors(int*& deviceFactorsCount, uint64_t *& deviceFactors, uint32_t deviceFactorsSize)
{
    cudaMalloc(&deviceFactorsCount, sizeof(*deviceFactorsCount));
    cudaMemset(deviceFactorsCount, 0, sizeof(*deviceFactorsCount));
    cudaMalloc(&deviceFactors, sizeof(*deviceFactors)*deviceFactorsSize);
}

void free_device_factors(int*& deviceFactorsCount, uint64_t *& deviceFactors)
{
    cudaFree(deviceFactorsCount);
    cudaFree(deviceFactors);
}


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
    out << "Usage: cutrialdive -i <index> --tf-bits <exponent> [--prp-test]" << std::endl;
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
    
    std::cout << "Using up to " << omp_get_max_threads() << " thread(s) on the CPU" << std::endl;

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
        int * deviceFactorsCount = nullptr;
        uint64_t * deviceFactors = nullptr;
        uint32_t const deviceFactorsSize = 65536;
        alloc_device_factors(deviceFactorsCount, deviceFactors, deviceFactorsSize);



        int numSMs;
        cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, 0);
        auto tfStart = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> sieve_time, computeDataTime, copyPrimeDataToHostTime;
    
        std::vector<uint64_t> primes;
        primes.reserve(1ull << 22);
        prime_data<uint64_t> primeData;
        primeData.reserve(1ull << 22);
        device_prime_data<uint64_t> devicePrimeData;

        for (auto k0 = uint64_t{2}, k1 = (std::min)(uint64_t{1} << 25, sieveUpperBound)
                ; k1 <= sieveUpperBound
                ; k0 = k1, k1 += uint64_t{1} << 25
            ) {
            auto sieveStart = std::chrono::high_resolution_clock::now();
            primes.clear();
            cutrialdive::sieve(k0, k1, primes);
            sieve_time += std::chrono::high_resolution_clock::now() - sieveStart;

            auto computeDataStart = std::chrono::high_resolution_clock::now();
            compute_prime_data(primes, primeData);
            computeDataTime += std::chrono::high_resolution_clock::now() - computeDataStart;

            device_synchronize();

            auto copyPrimeDataToHostStart = std::chrono::high_resolution_clock::now();
            devicePrimeData.copy_from_host(primes, primeData);
            copyPrimeDataToHostTime += std::chrono::high_resolution_clock::now() - copyPrimeDataToHostStart;
            
            if (k1 <= 2642258) {
                trial_div<3, uint64_t><<<256*numSMs, 256>>>(cn, nlen, devicePrimeData.get_data(), deviceFactorsCount, deviceFactors, deviceFactorsSize);
            } else if (k1 <= (1ull << 32)) {
                trial_div<2, uint64_t><<<256*numSMs, 256>>>(cn, nlen, devicePrimeData.get_data(), deviceFactorsCount, deviceFactors, deviceFactorsSize);
            } else {
                trial_div<1, uint64_t><<<256*numSMs, 256>>>(cn, nlen, devicePrimeData.get_data(), deviceFactorsCount, deviceFactors, deviceFactorsSize);
            }
        }
        device_synchronize();

        auto tfEnd = std::chrono::high_resolution_clock::now();
        

        std::cout << "Computing prime data on host took " << computeDataTime.count() << "s (cumulated time)" << std::endl;
        std::cout << "Copying prime data to device took " << copyPrimeDataToHostTime.count() << "s (cumulated time)" << std::endl;

        auto foundFactors = getFactors(sm671, deviceFactorsCount, deviceFactors, deviceFactorsSize);
        free_device_factors(deviceFactorsCount, deviceFactors);
        std::cout << "Found " << foundFactors.size() << " factor(s)";
        bool isFirst = true;
        for(auto const & f : foundFactors) {
            if(!isFirst) {
                std::cout << ", ";
            } else {
                std::cout << ":" << std::endl;
            }
            std::cout << f.first;
            if (f.second > 1) {
                std::cout << "^" << f.second;
            }
            isFirst = false;
        }

        std::cout << std::endl << "[Factoring took "
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

