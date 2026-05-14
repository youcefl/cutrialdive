/*
 * Creation date: 2026.04.21
 * Created by Youcef Lemsafer
 */
#pragma once

#include <cstdint>
#include <ostream>
#include "common_defines.h"
#include "barrett_mu_types.hpp"
#include "number_sequence.hpp"
#include "hgint.hpp"
#include "modular_arithmetic_detail.hpp"

namespace cutrialdive
{
    // Smarandache base 10 i.e. concatenation of the first n integers > 0 written in base 10.
    // S(1) = 1, S(2) = 12, ..., S(10) = 12345678910, ...
    template <>
    struct number_sequence<mode_flag::smarandache>
    {
        using index_type = uint64_t;
        using value_type = HgInt;
        using residue_type = uint64_t;
        using barrett_mu_type = barrett_mu_both_t;

        static char const* short_name();
        static value_type value(index_type n);
        static void print_value(index_type n, std::ostream & out);
        static void print_expression(index_type n, std::ostream & out);
        CUTRIALDIVE_DEVICE_AND_HOST static residue_type value_mod_2(index_type n);
        CUTRIALDIVE_DEVICE_AND_HOST static residue_type next_value_mod(residue_type v_n_mod_p, index_type n, residue_type p);
        CUTRIALDIVE_DEVICE_AND_HOST static residue_type next_value_mod_mu(residue_type v_n_mod_p, index_type n, residue_type p, barrett_mu_type mu_p);
        CUTRIALDIVE_DEVICE_AND_HOST static residue_type next_value_mod_2(residue_type v_n_mod_p, index_type n);
        /// Computes the shift by which S(n) needs to be multiplied in order to append n+1
        /// e.g. S(9) * shift + 10 = S(10), shift is 100, so the returned value is simply
        /// ceil(log10(n+1)).
        CUTRIALDIVE_DEVICE_AND_HOST static index_type compute_shift(index_type n);
    };
}


namespace cutrialdive {


    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    typename number_sequence<mode_flag::smarandache>::index_type
    number_sequence<mode_flag::smarandache>::compute_shift(index_type n)
    {
        auto m = n + 1;
        // Most frequent cases first
        if(m >= 10'000 && m < 100'000) {
            return 100'000;
        } else if(m >= 1000 && m < 10'000) {
            return 10'000;
        } else if(m >= 100'000 && m < 1'000'000) {
            return 1'000'000;
        } else if(m >= 1'000'000 && m < 10'000'000) {
            return 10'000'000;
        }
        // Fallback for all other cases
        index_type result = 10;
        while(m /= 10) {
            result *= 10;
        }
        return result;
    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    uint64_t number_sequence<mode_flag::smarandache>::value_mod_2(
        uint64_t n
    )
    {
        return n & 1;
    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    uint64_t number_sequence<mode_flag::smarandache>::next_value_mod(
        uint64_t v_n_mod_p,
        uint64_t n,
        uint64_t p
    )
    {
        return (__uint128_t(v_n_mod_p) * compute_shift(n) + n + 1) % p;
    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    uint64_t number_sequence<mode_flag::smarandache>::next_value_mod_mu(
        uint64_t v_n_mod_p,
        uint64_t n,
        uint64_t p,
        barrett_mu_type mu_p
    )
    {
        auto shift = compute_shift(n);
        auto a = __uint128_t(v_n_mod_p) * shift + n + 1;
        auto q = (a >> 64) ? mul128hi(a, mu_p.mu128) : ((a * mu_p.mu64) >> 64);
        __uint128_t r = a - q * p;
        return r >= p ? r - p : r;
    }
    
    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    uint64_t number_sequence<mode_flag::smarandache>::next_value_mod_2(
        uint64_t v_n_mod_p,
        uint64_t n
    )
    {
        return (n & 1) ^ 1;
    }


}
