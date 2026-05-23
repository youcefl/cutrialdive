/*
* Created on 2026.04.21
* Copyright (c) Youcef Lemsafer
*/
#pragma once

#include <cstdint>
#include <bit>

#include "common_defines.h"
#include "barrett_mu_types.hpp"
#include "modular_arithmetic_detail.hpp"
#include "hgint.hpp"


namespace cutrialdive {

    /// Smarandache base b, Sm_b(n) is the concatenation of the first n integers written in base b
    template <uint64_t Base>
    class smarandache
    {
        static_assert(Base > 1, "Base must be at least 2");
    public:
        using index_type = uint64_t;
        using value_type = HgInt;
        using residue_type = uint64_t;
        using barrett_mu_type = barrett_mu_both_t;

        static constexpr auto base  = Base;

        /// Returns a short name for the number sequence e.g. Sm3 for Smarandache base 3
        static char const* short_name();
        /// Returns Sm_b(n)
        static value_type value(index_type n);
        /// Prints Sm_b(n) to the given stream
        static void print_value(index_type n, std::ostream & out);
        /// Prints to the given stream an expression that evaluates to Sm_b(n)
        static void print_expression(index_type n, std::ostream & out);
        /// Returns Sm_b(n) mod 2
        CUTRIALDIVE_DEVICE_AND_HOST static residue_type value_mod_2(index_type n);
        /// Returns Sm_b(n + 1) mod p given Sm_b(n) mod p, n and p
        /// @pre p is odd
        CUTRIALDIVE_DEVICE_AND_HOST static residue_type next_value_mod(
            residue_type v_n_mod_p,
            index_type n,
            residue_type p
        );
        /// Returns Sm_b(n + 1) mod p given Sm_b(n) mod p, n, p and the Barrett reciprocal of p
        /// @pre p is odd
        CUTRIALDIVE_DEVICE_AND_HOST static residue_type next_value_mod_mu(
            residue_type v_n_mod_p,
            index_type n,
            residue_type p,
            barrett_mu_type mu_p
        );
        /// Returns Sm_b(n + 1) mod 2 given Sm_b(n) mod 2 and n
        CUTRIALDIVE_DEVICE_AND_HOST static residue_type next_value_mod_2(
            residue_type v_n_mod_p,
            index_type n
        );
    private:
        CUTRIALDIVE_DEVICE_AND_HOST static index_type compute_multiplier(index_type n);
        /// floor(log2(base))
        static constexpr int floor_log2_base = std::bit_width(base) - 1;
    };

}


// Impl
namespace cutrialdive {

    template <uint64_t Base>
    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    uint64_t smarandache<Base>::value_mod_2(
        uint64_t n
    )
    {
        if constexpr(!(base & 1)) {
            return n & 1;
        } else {
            // Sm_b(n) mod 2 = 1 if and only if n mod 4 belongs to {1, 2}
            return (n & 3) && ((n & 3) != 3);
        }
    }
    
    template <uint64_t Base>
    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    uint64_t smarandache<Base>::next_value_mod_2(
        uint64_t v_n_mod_p,
        uint64_t n
    )
    {
        if constexpr(!(base & 1)) {
            // For even bases the multiplier B^(logB(n + 1) + 1) is always even
            return (n & 1) ^ 1;
        } else {
            // v_n_mod_p is guaranteed to be a residue mod 2
            // the shift is a power of odd base thus always odd
            return v_n_mod_p ^ (n & 1) ^ 1;
        }
    }

    template <uint64_t Base>
    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    uint64_t smarandache<Base>::next_value_mod(
        uint64_t v_n_mod_p,
        uint64_t n,
        uint64_t p
    )
    {
        if constexpr(!(base & (base - 1))) {
            // When base is a power of 2
            auto shift = digits_in_base<base>(n + 1);
            return ((__uint128_t(v_n_mod_p) << (shift * floor_log2_base)) + n + 1) % p;
        } else {
            return (__uint128_t(v_n_mod_p) * compute_multiplier(n) + n + 1) % p;
        }
    }

    template <uint64_t Base>
    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    uint64_t smarandache<Base>::next_value_mod_mu(
        uint64_t v_n_mod_p,
        uint64_t n,
        uint64_t p,
        barrett_mu_type mu_p
    )
    {
        if constexpr(!(base & (base - 1))) {
            // When base is a power of 2
            auto shift = digits_in_base<base>(n + 1);
            auto a = (__uint128_t(v_n_mod_p) << (shift * floor_log2_base)) + n + 1;
            if(a >> 64) {
                auto q = mul128hi(a, mu_p.mu128);
                auto r = a - q * p;
                return r >= p ? r - p : r;
            }
            auto q = (__uint128_t(a) * mu_p.mu64) >> 64;
            auto r = a - q * p;
            return r >= p ? r - p : r;
        } else {
            auto a = __uint128_t(v_n_mod_p) * compute_multiplier(n) + n + 1;
            auto q = (a >> 64) ? mul128hi(a, mu_p.mu128) : ((a * mu_p.mu64) >> 64);
            __uint128_t r = a - q * p;
            return r >= p ? r - p : r;
        }
    }

    template <uint64_t Base>
    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    typename smarandache<Base>::index_type
    smarandache<Base>::compute_multiplier(index_type n)
    {
        // Note: this path is not taken when base is a power of two
        // because in that case we use bit shifts
        auto m = n + 1;
        // For some base we test the most frequent cases first
        if constexpr (base == 3) {
            if ((m >= 19683) && (m < 59049)) {
                return 59049;
            } else if((m >= 6561) && (m < 19683)) {
                return 19683;
            } else if ((m >= 177147) && (m < 531441)) {
                return 531441;
            } else if ((m >= 531441) && (m < 1594323)) {
                return 1594323;
            } else if ((m >= 1594323) && (m < 4782969)) {
                return 4782969;
            }
        } else if constexpr (base == 10) {
            if(m >= 10'000 && m < 100'000) {
                return 100'000;
            } else if(m >= 1000 && m < 10'000) {
                return 10'000;
            } else if(m >= 100'000 && m < 1'000'000) {
                return 1'000'000;
            } else if(m >= 1'000'000 && m < 10'000'000) {
                return 10'000'000;
            }
        }
        // Fallback for all other cases
        index_type result = base;
        while(m /= base) {
            result *= base;
        }
        return result;
    }
}
