/*
* Created on 2026.04.30
* Copyright Youcef Lemsafer
*/
#pragma once

#include <vector>
#include <cstdint>

#include "common_defines.h"

namespace cutrialdive {

    // Returns normalized d i.e. d*2^k such that 2^63 <= d*2^k < 2^64
    template <typename NumT>
    CUTRIALDIVE_DEVICE_AND_HOST NumT normalize(NumT d);

    // Returns floor((B^2 - 1) / d) - B = floor(((B - 1 - d) * B + B - 1) / d) where B = 2^64
    // Assumes 2^63 <= d < 2^64 i.e. most significant bit of d must be set
    template <typename NumT>
    CUTRIALDIVE_DEVICE_AND_HOST NumT reciprocal(NumT d);

    /// @brief Computes the reciprocal of each prime in @param primes and puts the result in @param reciprocals
    /// @param[in] primes vector of prime numbers
    /// @param[out] reciprocals vector such that reciprocals[i] == reciprocal(primes[i]) for each i < primes.size()
    void compute_reciprocals(std::vector<uint64_t> const & primes, std::vector<uint64_t> & reciprocals);

//#if CUTRIALDIVE_IS_CUDA

    template <typename NumT>
    struct initial_residue {
        NumT value;
        explicit CUTRIALDIVE_DEVICE_AND_HOST initial_residue(NumT);
    };

    // Computes N % d where N = <n[nlen-1], ..., n[1], n[0]>
    //  = n[0] + 2^64*n[1] + ... + 2^(64*(nlen-1))*n[nlen-1]
    CUTRIALDIVE_DEVICE uint64_t modnby1(uint64_t * n, uint64_t nlen, uint64_t d, initial_residue<uint64_t> r = initial_residue<uint64_t>{0});
    CUTRIALDIVE_DEVICE uint64_t modnby1(uint64_t * n, uint64_t nlen, uint64_t d, uint64_t rnd, initial_residue<uint64_t> r = initial_residue<uint64_t>{0});
    template <std::size_t N>
    CUTRIALDIVE_DEVICE uint64_t modnby1(uint64_t (&n) [N], uint64_t d, initial_residue<uint64_t> r = initial_residue<uint64_t>{0});

    // Returns (u1 * 2^64 + u0) mod d
    // Assumes u1 < d, 2^63 <= d < 2^64 and rnd is the reciprocal of d
    // i.e. rnd = floor((2^128 - 1) / d) - d
    CUTRIALDIVE_DEVICE uint64_t mod2by1(uint64_t u1, uint64_t u0, uint64_t d, uint64_t rnd);

//#endif // CUTRIALDIVE_IS_CUDA

    template <typename NumT>
    inline
    CUTRIALDIVE_DEVICE_AND_HOST
    initial_residue<NumT>::initial_residue(NumT val)
        : value(val)
    {
    }

    template <typename NumT>
    inline
    CUTRIALDIVE_DEVICE_AND_HOST
    NumT normalize(NumT d)
    {
#if CUTRIALDIVE_IS_CUDA
        // CUDA-specific __nv_clzll returns number of consecutive leading zero bits,
        // starting at most significant bit.
        return d << __clzll(d);
#else
        return d << __builtin_clzll(d);
#endif
    }

#if CUTRIALDIVE_IS_CUDA
    // Used inside reciprocal(uint64_t) (and only there).
    inline CUTRIALDIVE_DEVICE
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
#endif // CUTRIALDIVE_IS_CUDA

    // Returns floor((B^2 - 1) / d) - B = floor(((B - 1 - d) * B + B - 1) / d) where B = 2^64
    // Assumes 2^63 <= d < 2^64 i.e. most significant bit of d must be set
    template <typename NumT>
    inline
    CUTRIALDIVE_DEVICE_AND_HOST
    NumT reciprocal(NumT d)
    {
#if CUTRIALDIVE_IS_CUDA
        // <B-1-d, B-1> / d
        auto a = ~uint64_t{} - d;
        uint32_t a3 = a >> 32, a2 = a & 0xffffffff;
        uint32_t a1 = 0xffffffff, a0 = a1;
        auto q1 = inner_reciprocal(a3, a2, a1, d);
        auto q0 = inner_reciprocal(a2, a1, a0, d);
        return (uint64_t(q1) << 32) | q0;
#else
        static_assert(std::is_same_v<NumT, uint64_t>, "Integer type other than uint64_t not supported yet");
        if constexpr (std::is_same_v<NumT, uint64_t>) {
            constexpr auto all_ones = ~uint64_t{};
            return uint64_t(((__uint128_t(all_ones - d) << 64) | all_ones) / d);
        }
#endif
    }

#if CUTRIALDIVE_IS_CUDA

    // Returns (u1 * 2^64 + u0) mod d
    // Assumes u1 < d, 2^63 <= d < 2^64 and rnd is the reciprocal of d
    // i.e. rnd = floor((2^128 - 1) / d) - d
    inline
    CUTRIALDIVE_DEVICE uint64_t mod2by1(uint64_t u1, uint64_t u0, uint64_t d, uint64_t rnd)
    {
        uint64_t q1 = __umul64hi(u1, rnd), q0 = u1 * rnd;
        q0 += u0;
#if 0 // branchless version is not faster on GPU
        q1 += q0 < u0;
#else
        if (q0 < u0) {
            ++q1;
        }
#endif
        q1 += u1;
        ++q1;
        auto r = u0 - q1 * d;
#if 0 // branchless version is not faster on GPU
        r += uint64_t(-(r > q0)) & d;
        r -= uint64_t(-(r >= d)) & d;
#else        
        if (r > q0) {
            r += d;
        }
        if (r >= d) {
            r -= d;
        }
#endif
        return r;
    }

    inline
    CUTRIALDIVE_DEVICE uint64_t modnby1_internal(uint64_t * n, uint64_t nlen, uint64_t d, uint64_t nd, uint64_t rnd, uint64_t r)
    {
        for(auto j = nlen - 1; ; --j) {
            r = mod2by1(r, n[j], nd, rnd);
            if(!j) {
                break;
            }
        }
        return r % d;
    }

    // Computes N % d where N = <n[nlen-1], ..., n[1], n[0]>
    //  = n[0] + 2^64*n[1] + ... + 2^(64*(nlen-1))*n[nlen-1]
    inline
    CUTRIALDIVE_DEVICE uint64_t modnby1(uint64_t * n, uint64_t nlen, uint64_t d, initial_residue<uint64_t> ir)
    {
        auto nd = normalize(d);
        return modnby1_internal(n, nlen, d, nd, reciprocal(nd), ir.value);
    }

    inline
    CUTRIALDIVE_DEVICE uint64_t modnby1(uint64_t * n, uint64_t nlen, uint64_t d, uint64_t rnd, initial_residue<uint64_t> ir)
    {
        return modnby1_internal(n, nlen, d, normalize(d), rnd, ir.value);
    }

    template <size_t N>
    inline
    CUTRIALDIVE_DEVICE uint64_t modnby1(uint64_t (&n) [N], uint64_t d, initial_residue<uint64_t> ir)
    {
        return modnby1(&n[0], N, d, ir);
    }

#endif // CUTRIALDIVE_IS_CUDA

}
