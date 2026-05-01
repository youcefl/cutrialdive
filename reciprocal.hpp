/*
* Created on 2026.04.30
* Copyright Youcef Lemsafer
*/
#pragma once

#include <vector>
#include <cstdint>

namespace cutrialdive {

    // Returns normalized d i.e. d*2^k such that 2^63 <= d*2^k < 2^64
    template <typename NumT>
    inline
#ifdef __CUDA_ARCH__
    __device__ __host__
#endif
    NumT normalize(NumT d)
    {
    #if !defined __CUDA_ARCH__
        return d << __builtin_clzll(d);
    #else
        // CUDA-specific __nv_clzll returns number of consecutive leading zero bits,
        // starting at most significant bit.
        return d << __clzll(d);
    #endif
    }

#if defined __CUDA_ARCH__
    // Used inside reciprocal(uint64_t) (and only there).
    inline __device__
    uint_fast32_t inner_reciprocal(uint32_t & a2, uint32_t & a1, uint32_t & a0, uint64_t d)
    {
        auto q = ((uint64_t(a2) << 32) | a1) / (d >> 32);
        if (q > 0xffffffff) {
            q = 0xffffffff;
        }
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
#endif

    // Returns floor((B^2 - 1) / d) - B = floor(((B - 1 - d) * B + B - 1) / d) where B = 2^64
    // Assumes 2^63 <= d < 2^64 i.e. most significant bit of d must be set
    template <typename NumT>
#ifdef __CUDA_ARCH__
    __device__ __host__
#endif
    inline NumT reciprocal(NumT d)
    {
    #if !defined __CUDA_ARCH__
        if constexpr (std::is_same_v<NumT, uint64_t>) {
            constexpr auto all_ones = ~uint64_t{};
            return uint64_t(((__uint128_t(all_ones - d) << 64) | all_ones) / d);
        }
    #else  
        // <B-1-d, B-1> / d
        auto a = ~uint64_t{} - d;
        uint32_t a3 = a >> 32, a2 = a & 0xffffffff;
        uint32_t a1 = 0xffffffff, a0 = a1;
        auto q1 = inner_reciprocal(a3, a2, a1, d);
        auto q0 = inner_reciprocal(a2, a1, a0, d);
        return (uint64_t(q1) << 32) | q0;
    #endif
    }

    void compute_reciprocals(std::vector<uint64_t> const & primes, std::vector<uint64_t> & reciprocals);

}
