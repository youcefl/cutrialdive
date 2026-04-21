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

#include "hgint.hpp"

template <typename It>
std::ostream & output_range(std::ostream & out, It is, It ie)
{
    out << '{';
    auto sep = "";
    for(; is != ie; ++is) {
        out << sep << *is;
        sep = ", ";
    }
    out << '}';
    return out;
}

template <typename T>
std::ostream & output_range(std::ostream & out, T const & vec)
{
    return output_range(out, std::begin(vec), std::end(vec));
}



template <typename T>
auto smarandache(uint64_t n)
{
    if(!n || n == 1) {
        return T{n};
    }
    auto log10n = 0;
    for(uint64_t d = 1; (d *= 10) <= n;) {
        ++log10n;
    }
    const auto c3 = T{1490840987654358ull} * pow(T{10}, 180) - T{10990220};
    const auto c4 = pow(T{10}, 2701) * c3 * T{10201} - T{1109890222000ull};
    const auto c5 = pow(T{10}, 36000) * c4 * T{123454321} - T{123443211358} * pow(T{33300}, 2);
    static const std::array<std::function<T(uint64_t)>, 5> sm = {
          [](auto x) { return (pow(T{10}, x + 1) - T{9 * x + 10}) / T{81}; }
        , [](auto x) { return (pow(T{10}, 2*x - 18) * T{uint64_t(123456789)*99*99 + 991} 
                                                - T{99 * x + 100}) / T{9801}; }
        , [&](auto x) { return (c3 * pow(T{10}, 3*x - 296) - T{120879* x + 121000}) 
                                        / T{120758121}; }
        , [&](auto x) { return (c4 * pow(T{10}, 4*(x - 999)) - T{12321} * T{ 9999*x + 10000 }) 
                                        / T{1231853592321ull}; }
        , [&](auto x) { return (c5 * pow(T{10}, 5*(x - 9999)) - T{12321} * T{1234321} * (T{99999} * T{x} + T{100000})) 
                                        / (T{12321} * T{1234321} * T{99999} * T{99999}); }
        };
    if(log10n < sm.size()) {
        return sm[log10n](n);
    } else {
        throw std::exception{};
    }
}


template <typename It, typename Func, typename StopCond>
void parallel_for_each(int nbThreads, It first, It last, Func f, StopCond s)
{
    if (nbThreads < 2) {
        for (; first != last; ++first) {
            if (s(*first)) {
                break;
            }
            f(*first);
        }
        return;
    }
    std::vector<std::thread> workers;
    for(int i = 0; i < nbThreads; ++i) {
        if (std::distance(first, last) <= i) {
            break;
        }
        workers.emplace_back(std::thread{[&f, &s, first, last, i]() {
                const auto st = i + 1;
                for (auto it = first + i; ; it = it + st) {
                    if (s(*it)) {
                        break;
                    }
                    f(*it);
                    if (std::distance(it, last) <= st) {
                        break;
                    }
                } 
            }});
    }
    std::for_each(std::begin(workers), std::end(workers), [](auto & t) {
            if (t.joinable()) {
                t.join();
            }
        });
}

static const size_t DEVICE_FACTORS_SIZE = 20000;
__device__ uint64_t deviceFactors[DEVICE_FACTORS_SIZE];
__device__ int deviceFactorsCount = 0;


// Returns normalized d i.e. d*2^k such that 2^63 <= d*2^k < 2^64
__device__
uint64_t normalize(uint64_t d)
{
    // CUDA-specific __nv_clzll returns number of consecutive leading zero bits,
    // starting at most significant bit.
    return d << __clzll(d);
}


// Returns (u1 * 2^64 + u0) mod d
// Assumes u1 < d, 2^63 <= d < 2^64 and rnd is the reciprocal of d
// i.e. rnd = floor((2^128 - 1) / d) - d
__device__
uint64_t mod2by1(uint64_t u1, uint64_t u0, uint64_t d, uint64_t rnd)
{
    uint64_t q1 = __umul64hi(u1, rnd), q0 = u1 * rnd;
    q0 += u0;
    if (q0 < u0) {
        ++q1;
    }
    q1 += u1;
    ++q1;
    auto r = u0 - q1 * d;
    if (r > q0) {
        r += d;
    }
    if (r >= d) {
        r -= d;
    }
    return r;
}


// Used inside reciprocal(uint64_t) (and only there).
__device__
uint_fast32_t inner_reciprocal(uint32_t & a2, uint32_t & a1, uint32_t & a0, uint64_t d)
{
    auto q = ((uint64_t(a2) << 32) | a1) / (d >> 32);
    if (q > 0xffffffff)
        q = 0xffffffff;
    uint64_t qdl = q * d;
    uint64_t qdh = __umul64hi(q, d);
    while ((qdh > a2) || ((qdh == a2) && (qdl > ((uint64_t(a1) << 32) | a0))))  {
        --q;
        if (qdl < d) {
            --qdh;
        }
        qdl -= d;
    }
    if (a0 < (qdl & 0xffffffff)) {
        --a1;
    }
    a0 -= qdl & 0xffffffff;
    if (a1 < (qdl >> 32)) {
        --a2;
    }
    a1 -= qdl >> 32;
    a2 -= qdh;
    return q;        
}


// Returns floor((B^2 - 1) / d) - B = floor(((B - 1 - d) * B + B - 1) / d) where B = 2^64
// Assumes 2^63 <= d < 2^64 i.e. most significant bit of d must be set
__device__
uint64_t reciprocal(uint64_t d)
{
    // <B-1-d, B-1> / d
    auto a = ~uint64_t{} - d;
    uint32_t a3 = a >> 32, a2 = a & 0xffffffff;
    uint32_t a1 = 0xffffffff, a0 = a1;
    auto q1 = inner_reciprocal(a3, a2, a1, d);
    auto q0 = inner_reciprocal(a2, a1, a0, d);
    return (uint64_t(q1) << 32) | q0;
}

// Computes N % d where N = <n[nlen-1], ..., n[1], n[0]>
//  = n[0] + 2^64*n[1] + ... + 2^(64*(nlen-1))*n[nlen-1]
__device__
uint64_t modnby1(uint64_t * n, uint64_t nlen, uint64_t d)
{
    uint64_t r = 0;
    auto nd = normalize(d);
    auto rnd = reciprocal(normalize(d));
    for(auto j = nlen - 1; ; --j) {
        r = mod2by1(r, n[j], nd, rnd);
        if(!j) {
            break;
        }
    }
    return r % d;
}

template <size_t N>
__device__
uint64_t modnby1(uint64_t (&n) [N], uint64_t d)
{
    return modnby1(&n[0], N, d);
}


__device__
void pushFactor(uint64_t ff)
{
    auto factorIdx = atomicAdd(&deviceFactorsCount, 1);
    deviceFactors[factorIdx] = ff;
}

template <uint8_t BatchSize>
__global__
void trial_div(uint64_t * n, uint64_t nlen, uint64_t * p, size_t plen)
{
    int index = (threadIdx.x + blockIdx.x * blockDim.x) * BatchSize;
    int stride = blockDim.x * gridDim.x * BatchSize;

    for(size_t i = index; i < plen; i += stride) {
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
    printf("Adding factor %llu", p[i+j]);
                pushFactor(p[i + j]);
            }
        }
    }
}

std::vector<uint64_t> getFactorsFromDevice()
{
    //cudaMemcpy
    int factorsCount = 0;
    cudaMemcpy(&factorsCount, &deviceFactorsCount, sizeof(int), cudaMemcpyDeviceToHost);
    std::cout << "factorsCount = " << factorsCount << std::endl;
    auto ret = std::vector<uint64_t>(std::size_t(factorsCount), uint64_t(0));
    std::cout << "Size of ret = " << ret.size() << std::endl;
    cudaMemcpy(&ret[0], deviceFactors, sizeof(uint64_t)*factorsCount, cudaMemcpyDeviceToHost);
    std::sort(std::begin(ret), std::end(ret));
    return ret;
}


// Tests reciprocal on the device
__global__
void treciprocal()
{
    struct testcase { uint64_t input; uint64_t expected; };
    constexpr auto testcaseSize = 24;
    testcase tvec[testcaseSize] = { 
          {0x8000000000000000, 0xffffffffffffffff}, {0x8000000000000001, 0xfffffffffffffffc}
        , {0x8000000000000002, 0xfffffffffffffff8}, {0x8000000000000003, 0xfffffffffffffff4}
        , {0xffffffffffffffff, 0x0000000000000001}, {0xfffffffffffffffe, 0x0000000000000002}
        , {0xfffffffffffffffd, 0x0000000000000003}, {0xfffffffffffffffc, 0x0000000000000004}
        , {0xfffffffffffffffb, 0x0000000000000005}, {0xff00ff00ff00ff00, 0x0100000000000001}
        , {0xfe000000ffffffff, 0x0204080f1c3460b3}, {0xff00000000000000, 0x0101010101010101}
        , {0xfdfdfdfdfdfdfdfd, 0x02061236a3ebc34a}, {0xa68d164e738685dc, 0x897d219bfdaf4388}
        , {0x82bf6e750f72fbeb, 0xf53d68328e2b07e2}, {0x915be40a6e3459d6, 0xc2db378ada15a0be}
        , {0xa890630f91cbdf17, 0x84ca3ca49f815c02}, {0x9603e01e9352f3c9, 0xb4dcd1d74a60f12a}
        , {0xb88072660b8881aa, 0x63348f4a3d64d9f7}, {0xfbbc44125f5cb984, 0x04563a985b4a878f}
        , {0x935613859937496c, 0xbcce438cd4fbc361}, {0xf1797c593356ac16, 0x0f6631c930b05025}
        , {0xe28323756b21a609, 0x215395b2245b5ddd}, {0x9541d0d870c791c3, 0xb714d0c552f0a5c9}
    };

    for (auto i = 0; i < testcaseSize; ++i) {
        auto const & t = tvec[i];
        auto expect = t.expected;
        auto actual = reciprocal(t.input);
        if (expect != actual) {
            printf("FAIL -> reciprocal(%" PRIu64 ")\n", t.input); 
            printf("    expected: %" PRIu64 " | <low = %u, high = %u>\n", expect, uint32_t(expect & 0xffffffff), uint32_t(expect >> 32));
            printf("      actual: %" PRIu64 " | <low = %u, high = %u>\n", actual, uint32_t(actual & 0xffffffff), uint32_t(actual >> 32));
        }
    }
}

// Test mod2by1 on the device
__global__
void tmod2by1()
{
    struct testcase { uint64_t u1, u0, d, rnd, expected; };
    testcase tvec[] = { 
          {               0x0,                0x0, 0x8000000000000000, 0xffffffffffffffff,                0x0}
        , {               0x1,                0x1, 0x8000000000000000, 0xffffffffffffffff,                0x1}
        , {0x0fffffffffffffff, 0x0fffffffffffffff, 0x8000000000000000, 0xffffffffffffffff, 0x0fffffffffffffff}
        , {               0x0, 0xffffffffffffffff, 0xc5bb533e4317e471, 0x4b705f0c618a3b59, 0x3a44acc1bce81b8e}
        , {0xfffffffffffffffe, 0xffffffffffffffff, 0xffffffffffffffff, 0x0000000000000001, 0xfffffffffffffffe}
        , {0xfcfcfcfcfcfcfcfc, 0xffffffffffffffff, 0xfcfcfcfcfcfcfcfd, 0x030c30c30c30c30c, 0xfcfcfcfcfcfcfcfc}
        , {        0xffffffff, 0xeeeeeeeeeeeeeeee, 0x97fe9508a59e5231, 0xaf2c716780c9f84d, 0x3b6623c2ca883c37}
    };
    constexpr auto testcaseSize = sizeof(tvec)/sizeof(tvec[0]);

    for (auto i = 0; i < testcaseSize; ++i) {
        auto const & t = tvec[i];
        auto expect = t.expected;
        auto actual = mod2by1(t.u1, t.u0, t.d, t.rnd);
        if (expect != actual) {
            printf("FAIL -> mod2by1(<%" PRIu64 ", %" PRIu64 ">, %" PRIu64 ", %" PRIu64 ")\n", t.u1, t.u0, t.d, t.rnd); 
            printf("    expected: %" PRIu64 " | <low = %u, high = %u>\n", expect, uint32_t(expect & 0xffffffff), uint32_t(expect >> 32));
            printf("      actual: %" PRIu64 " | <low = %u, high = %u>\n", actual, uint32_t(actual & 0xffffffff), uint32_t(actual >> 32));
        }
    }
}


__global__
void tmodnby1()
{
    uint64_t two192m1[] = {0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff};
    uint64_t m521[] = { 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff
                      , 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff
                      , 511 };

    struct testcase { uint64_t * n, nlen, d, expected; };
#define NS(n) n, sizeof(n)/sizeof(n[0])
#define TC testcase
    testcase tvec[] = {
          TC{NS(two192m1), 1, 0}, TC{NS(two192m1), 2, 1}, TC{NS(two192m1), 3, 0}
        , TC{NS(two192m1), 197, 36}, TC{NS(two192m1), 0xffffffffffffffff, 0}
        , TC{NS(m521), 2, 1}, TC{NS(m521), 1 << 13, (1 << 13) - 1}
        , TC{NS(m521), 1ull << 63, (1ull << 63) - 1}
    };
#undef NS
#undef TC  
    constexpr auto testcaseSize = sizeof(tvec)/sizeof(tvec[0]);
    
    for (auto i = 0; i < testcaseSize; ++i) {
        auto const & t = tvec[i];
        auto expect = t.expected;
        auto actual = modnby1(t.n, t.nlen, t.d);
        if (expect != actual) {
            printf("FAIL -> modnby1(2^192-1, 2)\n");
            printf("    expected: %" PRIu64 "\n", expect);
            printf("      actual: %" PRIu64 "\n", actual);
        }
    }
    
}


const std::array<uint64_t, 54> primes8 = {
      2,   3,   5,   7,  11,  13,  17,  19,  23,  29,
     31,  37,  41,  43,  47,  53,  59,  61,  67,  71,
     73,  79,  83,  89,  97, 101, 103, 107, 109, 113,
    127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
    179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
    233, 239, 241, 251
};


const std::vector<uint64_t> primes16 = [](){
        std::vector<uint64_t> primes;
        primes.reserve(6542); // pi(2^16) = 6542
        std::copy(std::begin(primes8), std::end(primes8), std::back_inserter(primes));
        uint8_t wheel[] = {2, 4};
        uint8_t i = 0;
        for (uint64_t n = 257; n < (1 << 16); n += wheel[i], i ^= 1) {
            auto isprime = true;
            for (auto j = 1; j < primes8.size(); ++j) {
                auto d = primes8[j];
                if (!(n % d)) {
                    isprime = false;
                    break;
                }
                if ( d * d > n ) {
                    break;
                }
            }
            if (isprime) {
                primes.push_back(n);
            }
        }
        return primes;
    }();

uint32_t bitwidth(uint64_t t)
{
    uint32_t c{};
    while (t >>= 1) {
        ++c;
    }
    return c;
}

uint64_t intsqrt(uint64_t n)
{
    if (n < 2) {
        return n;
    }
    auto x0 = uint64_t(1) << (bitwidth(n) / 2 + 1);
    for (auto x1 = (x0 + n / x0) / 2
        ; x1 < x0
        ; x0 = x1, x1 = (x0 + n / x0) / 2
        ) {
    }
    return x0;
}


void sieve(std::vector<uint64_t> & primes, uint64_t n0, uint64_t n1)
{
    assert(n0 <= n1);
    primes.clear();
    if(n1 <= (1 << 16)) {
        for(auto p : primes16) {
            if(n0 <= p) {
                if(p < n1) {
                    primes.push_back(p);
                } else {
                    break;
                }
            }
        }
        return;
    }
    if (n1 <= (uint64_t(1) << 32)) {
        if (n0 < (1 << 16)) {
            auto it = std::lower_bound(std::begin(primes16), std::end(primes16), n0);
            primes.reserve((n1 - (1 << 16) + 1) / 2 + std::distance(it, std::end(primes16)));
            std::copy(it, std::end(primes16), std::back_inserter(primes));
            n0 = 1 << 16;
        } else {
            primes.reserve((n1 - n0 + 1) / 2);
        }
        auto i0 = primes.size();
        n0 = n0 + ((n0 ^ 1) & 1);
        for (auto x = n0; x < n1; x += 2) {
            primes.push_back(x);
        }
        if (primes.size() == i0) {
            return;
        }
        /*parallel_for_each(4, std::begin(primes16) + 1, std::end(primes16), [&primes, n0, i0](uint64_t sp) {
                    auto s = ((n0 + sp - 1) / sp) * sp;
                    for(auto k = i0 + ((s & 1) ? (s - n0) / 2 : (s + sp - n0) / 2); k < primes.size(); k += sp) {
                        primes[k] = 0;
                    }
              }
            , [n1](uint64_t x) { return x * x >= n1; }
            );*/
        for (auto i = 1; i < primes16.size(); ++i) {
            auto sp = primes16[i];
            if( sp * sp >= n1 ) {
                break;
            }
            auto s = ((n0 + sp - 1) / sp) * sp;
            for(auto k = i0 + ((s & 1) ? (s - n0) / 2 : (s + sp - n0) / 2); k < primes.size(); k += sp) {
                primes[k] = 0;
            }
        }
        primes.erase(std::remove(std::begin(primes) + i0, std::end(primes), 0), std::end(primes));
        return;
    }
    if (n0 < uint64_t(1) << 32) {
        sieve(primes, n0, uint64_t(1) << 32);
        n0 = uint64_t(1) << 32;
    }
    auto sieves = primes.size();
    primes.reserve(sieves + (n1 - n0 + 1) / 2);
    n0 = n0 + ((n0 ^ 1) & 1);
    for (auto x = n0; x < n1; x += 2) {
        primes.push_back(x);
    }
    if (primes.size() == sieves) {
        return;
    }
    // Following loop can be parallelized, idea: construct vector of pairs v = [start, end)
    // then std::for_each(std::par, v, []() { sieve });
    // of course if doing so tmpPrimes can no longer be shared among iterations of the loop.
    constexpr auto sieveRangeLength = uint64_t(1) << 21;
    const auto upperBound = intsqrt(n1) + 1;
    std::vector<uint64_t> tmpPrimes;
    for ( auto ks = uint64_t(3)
        ; ks < upperBound
        ; ks += sieveRangeLength
        ) {
        auto ke = (std::min)(ks + sieveRangeLength, upperBound);
        sieve(tmpPrimes, ks, ke);
        for (auto p : tmpPrimes) {
            auto i0 = sieves;
            if (p * p >= n0) {
                i0 += (p * p - n0) / 2;
            } else {
                auto s = ((n0 + p - 1) / p) * p;
                i0 += (s & 1) ? (s - n0) / 2 : (s + p - n0) / 2;
            }
            for (auto i = i0; i < primes.size(); i += p) {
                primes[i] = 0;
            }
        }
    }
    primes.erase(std::remove(std::begin(primes) + sieves, std::end(primes), 0), std::end(primes));
}

std::vector<uint64_t> sieve(uint64_t n0, uint64_t n1)
{
    std::vector<uint64_t> primes;
    sieve(primes, n0, n1);
    return primes;    
}


int autotest()
{
        std::cout << "Executing self tests..." << std::endl;
        {
            std::vector<std::pair<uint64_t,uint64_t>> tsqrt = { {0, 0}, {1, 1}, {2, 1}, {3, 1}, {4, 2}
                , {17, 4}, {1024, 32}, {2047, 45}, {2048, 45}, {2115, 45}, {2116, 46}, {~uint64_t{}, 4294967295ull} };
            for (auto const & t : tsqrt) {
                auto sqrtt = intsqrt(t.first);
                if(sqrtt != t.second) {
                    std::cerr << "Failure: expected integer square root of " << t.first << " is " << t.second
                      << ", actual: " << sqrtt << std::endl;
                }
            }
        }
        {
            auto primes = sieve(256, 307);
            std::vector<uint64_t> expected = {257, 263, 269, 271, 277,
                                              281, 283, 293};
            if (primes != expected) {
                std::cerr << "Failure: expected 8 primes in [256, 307)" << std::endl;
                std::cerr << "Expected: "; output_range(std::cerr, expected); std::cerr << std::endl;
                std::cerr << "  Actual: "; output_range(std::cerr, primes); std::cerr << std::endl;
            }
        }
        {
            auto primes = sieve(1 << 16, (1 << 16) + 1);
            if (primes.size()) {
                std::cerr << "Failure: expected no primes in [65536, 65537)" << std::endl;
            }
        }
        {
            auto primes = sieve(1 << 16, (1 << 16) + 2);
            if (primes.size() != 1 || primes.front() != 65537) {
                std::cerr << "Failure: expected 65537 as the unique prime in [65536, 65538)" << std::endl;
                std::cerr << "size: " << primes.size() << std::endl;
            }
        }
        {
            auto primes = sieve(1 << 16, (1 << 16) + 100);
            std::vector<uint64_t> expected = {65537, 65539, 65543, 65551, 65557,
                                              65563, 65579, 65581, 65587, 65599,
                                              65609, 65617, 65629, 65633};
            if (primes != expected) {
                std::cerr << "Failure: expected 14 primes in [65536, 65536+100)" << std::endl;
                std::cerr << "Expected: "; output_range(std::cerr, expected); std::cerr << std::endl;
                std::cerr << "  Actual: "; output_range(std::cerr, primes); std::cerr << std::endl;
            }
        }
        {
            auto primes = sieve((1ull << 31) - 1, (1ull << 31) + 100);
            std::vector<uint64_t> expected = {2147483647, 2147483659, 2147483693,
                                              2147483713, 2147483743};
            if (primes != expected) {
                std::cerr << "Failure: expected 5 primes in [2^31-1, 2^31+100)" << std::endl;
                std::cerr << "Expected: "; output_range(std::cerr, expected); std::cerr << std::endl;
                std::cerr << "  Actual: "; output_range(std::cerr, primes); std::cerr << std::endl;
            }
        }
        {
            auto primes = sieve(3299100100, 3299100100 + 256);
            std::vector<uint64_t> expected = {3299100107, 3299100113, 3299100127,
                                              3299100191, 3299100199, 3299100221,
                                              3299100239, 3299100241, 3299100253,
                                              3299100269, 3299100293, 3299100317,
                                              3299100319, 3299100347};
            if (primes != expected) {
                std::cerr << "Failure: expected 14 primes in [3299100100, 3299100100 + 256)" << std::endl;
                std::cerr << "Expected: "; output_range(std::cerr, expected); std::cerr << std::endl;
                std::cerr << "  Actual: "; output_range(std::cerr, primes); std::cerr << std::endl;
            }
        }
        {
            auto primes = sieve(1, 1 << 17);
            const auto expectedCount = 12251;
            if (primes.size() != expectedCount) {
                std::cerr << "Failure: expected number of primes in [1, 131072)" << std::endl;
                std::cerr << "Expected: " << expectedCount << std::endl;
                std::cerr << "  Actual: " << primes.size() << std::endl;
            }
            if (!std::equal(std::begin(primes8), std::end(primes8), std::begin(primes))) {
                std::cerr << "Failure: " << std::endl;
                std::cerr << "Expected: "; output_range(std::cerr, primes8); std::cerr << std::endl;
                std::cerr << "  Actual: "; output_range(std::cerr, std::begin(primes)
                                            , std::begin(primes) + primes8.size()); std::cerr << std::endl;
            }
        }
        {
            auto primes = sieve(1ull << 32, (1ull << 32) + 200);
            const auto expectedCount = 8;
            if (primes.size() != expectedCount) {
                std::cerr << "Failure: expected number of primes in [2^32, 2^32+8)" << std::endl;
                std::cerr << "Expected: " << expectedCount << std::endl;
                std::cerr << "  Actual: " << primes.size() << std::endl;
            }
            std::vector<uint64_t> expected = {4294967311ull, 4294967357ull, 4294967371ull, 4294967377ull,
                                              4294967387ull, 4294967389ull, 4294967459ull, 4294967477ull};
            if (!std::equal(std::begin(primes), std::end(primes)
                , std::begin(expected), std::end(expected))) {
                std::cerr << "Failure: " << std::endl;
                std::cerr << "Expected: "; output_range(std::cerr, expected); std::cerr << std::endl;
                std::cerr << "  Actual: "; output_range(std::cerr, primes); std::cerr << std::endl;
            }
        }
        {
            auto primes = sieve(1ull << 32, (1ull << 32) + 1000);
            const auto expectedCount = 56;
            if (primes.size() != expectedCount) {
                std::cerr << "Failure: expected number of primes in [2^32, 2^32+1000)" << std::endl;
                std::cerr << "Expected: " << expectedCount << std::endl;
                std::cerr << "  Actual: " << primes.size() << std::endl;
            }
            std::vector<uint64_t> expected = {4294967311, 4294967357, 4294967371, 4294967377, 4294967387,
                                              4294967389, 4294967459, 4294967477, 4294967497, 4294967513,
                                              4294967539, 4294967543, 4294967549, 4294967561, 4294967563,
                                              4294967569, 4294967597, 4294967627, 4294967639, 4294967653,
                                              4294967681, 4294967687, 4294967701, 4294967723, 4294967759,
                                              4294967779, 4294967783, 4294967791, 4294967797, 4294967821,
                                              4294967857, 4294967861, 4294967867, 4294967891, 4294967969,
                                              4294967977, 4294967983, 4294968001, 4294968017, 4294968019,
                                              4294968023, 4294968059, 4294968101, 4294968127, 4294968143,
                                              4294968149, 4294968173, 4294968187, 4294968199, 4294968211,
                                              4294968233, 4294968239, 4294968257, 4294968283, 4294968287,
                                              4294968289};
            if (!std::equal(std::begin(primes), std::end(primes)
                , std::begin(expected), std::end(expected))) {
                std::cerr << "Failure: " << std::endl;
                std::cerr << "Expected: "; output_range(std::cerr, expected); std::cerr << std::endl;
                std::cerr << "  Actual: "; output_range(std::cerr, primes); std::cerr << std::endl;
            }
        }
        {
            auto primes = sieve((1ull << 32) - 256, (1ull << 32) + 256);
            const auto expectedCount = 22;
            if (primes.size() != expectedCount) {
                std::cerr << "Failure: expected number of primes in [2^32-256, 2^32+256)" << std::endl;
                std::cerr << "Expected: " << expectedCount << std::endl;
                std::cerr << "  Actual: " << primes.size() << std::endl;
            }
            std::vector<uint64_t> expected = {4294967087, 4294967111, 4294967143, 4294967161, 4294967189,
                                              4294967197, 4294967231, 4294967279, 4294967291, 4294967311,
                                              4294967357, 4294967371, 4294967377, 4294967387, 4294967389,
                                              4294967459, 4294967477, 4294967497, 4294967513, 4294967539,
                                              4294967543, 4294967549};
            if (!std::equal(std::begin(primes), std::end(primes)
                , std::begin(expected), std::end(expected))) {
                std::cerr << "Failure: " << std::endl;
                std::cerr << "Expected: "; output_range(std::cerr, expected); std::cerr << std::endl;
                std::cerr << "  Actual: "; output_range(std::cerr, primes); std::cerr << std::endl;
            }
        }
        {
            auto primes = sieve(65537ull*65537 - 200, 65537ull*65537 + 1);
            const auto expectedCount = 7;
            if (primes.size() != expectedCount) {
                std::cerr << "Failure: expected number of primes in [65537^2 - 200, 65537^2 + 1)" << std::endl;
                std::cerr << "Expected: " << expectedCount << std::endl;
                std::cerr << "  Actual: " << primes.size() << std::endl;
            }
            std::vector<uint64_t> expected = {4295098183, 4295098231, 4295098237, 4295098277, 4295098303,
                                              4295098309, 4295098349};
            if (!std::equal(std::begin(primes), std::end(primes)
                , std::begin(expected), std::end(expected))) {
                std::cerr << "Failure: " << std::endl;
                std::cerr << "Expected: "; output_range(std::cerr, expected); std::cerr << std::endl;
                std::cerr << "  Actual: "; output_range(std::cerr, primes); std::cerr << std::endl;
            }
        }
        {
            auto primes = sieve((1ull << 34) - 200, (1ull << 34) + 200);
            const auto expectedCount = 11;
            if (primes.size() != expectedCount) {
                std::cerr << "Failure: expected number of primes in [2^34-200, 2^34+200)" << std::endl;
                std::cerr << "Expected: " << expectedCount << std::endl;
                std::cerr << "  Actual: " << primes.size() << std::endl;
            }
            std::vector<uint64_t> expected = {17179868999, 17179869019, 17179869041, 17179869053,
                                              17179869071, 17179869107, 17179869143, 17179869209,
                                              17179869263, 17179869269, 17179869337};
            if (!std::equal(std::begin(primes), std::end(primes)
                , std::begin(expected), std::end(expected))) {
                std::cerr << "Failure: " << std::endl;
                std::cerr << "Expected: "; output_range(std::cerr, expected); std::cerr << std::endl;
                std::cerr << "  Actual: "; output_range(std::cerr, primes); std::cerr << std::endl;
            }
        }
        {
            auto s = (uint64_t(1) << 41)- (uint64_t(1) << 37) + (uint64_t(1) << 32);
            auto primes = sieve(s - 1, s + 1000);
            const auto expectedCount = 26;
            if (primes.size() != expectedCount) {
                std::cerr << "Failure: expected number of primes in [2^41-2^37+2^32-1, 2^41-2^37+2^32+1000)" << std::endl;
                std::cerr << "Expected: " << expectedCount << std::endl;
                std::cerr << "  Actual: " << primes.size() << std::endl;
            }
            std::vector<uint64_t> expected = {
                    2065879269491, 2065879269521, 2065879269527, 2065879269547, 2065879269563,
                    2065879269569, 2065879269617, 2065879269647, 2065879269749, 2065879269841,
                    2065879269869, 2065879269871, 2065879269899, 2065879269989, 2065879269991,
                    2065879270009, 2065879270057, 2065879270073, 2065879270093, 2065879270181,
                    2065879270201, 2065879270219, 2065879270249, 2065879270307, 2065879270327,
                    2065879270331};
            if (!std::equal(std::begin(primes), std::end(primes)
                , std::begin(expected), std::end(expected))) {
                std::cerr << "Failure: " << std::endl;
                std::cerr << "Expected: "; output_range(std::cerr, expected); std::cerr << std::endl;
                std::cerr << "  Actual: "; output_range(std::cerr, primes); std::cerr << std::endl;
            }
        }
        {
            auto primes = sieve((1ull << 61) - 100, (1ull << 61));
            const auto expectedCount = 3;
            if (primes.size() != expectedCount) {
                std::cerr << "Failure: expected number of primes in [2^61-100, 2^61-1)" << std::endl;
                std::cerr << "Expected: " << expectedCount << std::endl;
                std::cerr << "  Actual: " << primes.size() << std::endl;
            }
            std::vector<uint64_t> expected = {
                    2305843009213693907ull, 2305843009213693921ull, 2305843009213693951ull};
            if (!std::equal(std::begin(primes), std::end(primes)
                , std::begin(expected), std::end(expected))) {
                std::cerr << "Failure: " << std::endl;
                std::cerr << "Expected: "; output_range(std::cerr, expected); std::cerr << std::endl;
                std::cerr << "  Actual: "; output_range(std::cerr, primes); std::cerr << std::endl;
            }
        }
        {
            treciprocal<<<1, 1>>>();
            cudaDeviceSynchronize();
        }
        {
            tmod2by1<<<1, 1>>>();
            cudaDeviceSynchronize();
        }
        {
            tmodnby1<<<1, 1>>>();
            cudaDeviceSynchronize();
        }
        std::cout << "Done" << std::endl;
        return 0;
}


// Gets the factors from the device, then repeatedly divides the number by the factors
// to get the cofactor.
// number <- number / (f_1^e_1*f_2^e_2...f_n^e_n)
// @param [in/out] number which was sent to the GPU for factorization
// @return factors (with their exponent)
std::vector<std::pair<uint64_t, uint32_t>> getFactors(HgInt & number)
{
    std::vector<std::pair<uint64_t, uint32_t>> factors;
    auto const & rawFactors = getFactorsFromDevice();
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

extern bool isPRP(HgInt number);


int main(int argc, char** argv)
{
    if ((argc == 2) && (std::string(argv[1]) == "-t")) {
        return autotest();
    }
    uint64_t ii = 0;
    if ((argc == 3) && (std::string(argv[1]) == "-i")) {
        std::istringstream istr{argv[2]};
        istr >> ii;
    } else {
        std::cerr << "Usage: " << argv[0] << " -i <index>" << std::endl;
        return 1;
    }
    

    try {
        auto sm671 = smarandache<HgInt>(ii);
        auto * n = sm671.limbs();
        auto nlen = sm671.size();
        auto sizeInBase10 = sm671.sizeInBase(10);

        std::cout << "Smarandache(" << ii << ") = 0x" << std::hex << n[nlen-1] 
            << "..." << n[0] <<  ", length = " << std::dec << nlen << " 64-bit words, "
            <<  sizeInBase10 << " digits" << std::endl;

        uint64_t * cn;
        auto cuStatus = cudaMalloc(&cn, nlen*sizeof(uint64_t));
        if(cuStatus != cudaSuccess) {
            std::cerr << "Error returned by cudaMalloc(): " << cudaGetErrorString(cuStatus) << std::endl;
        }
        cudaMemcpy(cn, n, nlen*sizeof(uint64_t), cudaMemcpyHostToDevice);

        int numSMs;
        cudaDeviceGetAttribute(&numSMs, cudaDevAttrMultiProcessorCount, 0);
        auto tfStart = std::chrono::high_resolution_clock::now();
    
        std::vector<uint64_t> primes;
        uint64_t * devicePrimes{};
        auto primesSize = primes.size();
        for (auto k0 = uint64_t(2), k1 = uint64_t(1) << 21
                ; k1 <= (uint64_t(1) << 32)
                ; k0 = k1, k1 += uint64_t(1) << 21
            ) {
            sieve(primes, k0, k1);

            if (primes.size() > primesSize) {
                if (devicePrimes) {
                    cudaFree(devicePrimes);
                }
                cudaMalloc(&devicePrimes, primes.size() * sizeof(primes[0]));
            }
            primesSize = primes.size();
            cuStatus = cudaMemcpy(devicePrimes, &primes[0], primesSize * sizeof(primes[0]), cudaMemcpyHostToDevice);
            if (k1 <= 2642258) {
                trial_div<3><<<32*numSMs, 256>>>(cn, nlen, devicePrimes, primesSize);
            } else if (k1 <= (1ull << 32)) {
                trial_div<2><<<32*numSMs, 256>>>(cn, nlen, devicePrimes, primesSize);
            } else {
                trial_div<1><<<32*numSMs, 256>>>(cn, nlen, devicePrimes, primesSize);
            }
            cuStatus = cudaDeviceSynchronize();
            if(cuStatus != cudaSuccess) {
                std::cerr << "Error returned by cudaDeviceSynchronize(): " <<
                    cudaGetErrorString(cuStatus) << std::endl;
                int devCount = 0;
                cuStatus = cudaGetDeviceCount(&devCount);
                if(cuStatus != cudaSuccess) {
                    std::cerr << "Error returned by cudaGetDeviceCount(int*): " <<
                        cudaGetErrorString(cuStatus) << std::endl;
                }
            }
        }
        if (devicePrimes) {
            cudaFree(devicePrimes);
        }
        cudaFree(cn);
        auto tfEnd = std::chrono::high_resolution_clock::now();
        
        auto foundFactors = getFactors(sm671);
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

        std::cout << " [Factoring... "
                << std::chrono::duration<double, std::milli>(tfEnd - 
                        tfStart).count() / 1000 << "s]";
        std::flush(std::cout);
        auto prpStart = std::chrono::high_resolution_clock::now();
        auto hasPrpCofactor = isPRP(sm671);
        auto prpEnd = std::chrono::high_resolution_clock::now();
        std::cout << " [Primality testing... "
                << std::chrono::duration<double, std::milli>(prpEnd - 
                        prpStart).count() / 1000 << "s]"
                << std::endl;

        std::cout << "Smarandache(" << ii << ") -> " << (hasPrpCofactor ? "PRP" : "C") << std::endl;
        
    } catch (std::exception const & ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
        return 3;
    }
    return 0;
}

