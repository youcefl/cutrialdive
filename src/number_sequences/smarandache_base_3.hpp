/*
 * Creation date: 2026.05.19
 * Created by Youcef Lemsafer
 */
#pragma once

#include "common_defines.h"
#include "barrett_mu_types.hpp"
#include "number_sequence.hpp"
#include "hgint.hpp"
#include "modular_arithmetic_detail.hpp"

namespace cutrialdive {

    template <>
    struct number_sequence<num_seq_id::smarandache, 3>
    {
        static constexpr uint32_t base = 3;
        using index_type = uint64_t;
        using value_type = HgInt;
        using residue_type = uint64_t;
        using barrett_mu_type = barrett_mu_both_t;

        static char const* short_name();
        static value_type value(index_type n);
        static void print_value(index_type n, std::ostream & out);
        static void print_expression(index_type n, std::ostream & out);
        CUTRIALDIVE_DEVICE_AND_HOST static residue_type value_mod_2(index_type n);
        CUTRIALDIVE_DEVICE_AND_HOST static residue_type next_value_mod(
            residue_type v_n_mod_p,
            index_type n,
            residue_type p
        );
        CUTRIALDIVE_DEVICE_AND_HOST static residue_type next_value_mod_mu(
            residue_type v_n_mod_p,
            index_type n,
            residue_type p,
            barrett_mu_type mu_p
        );
        CUTRIALDIVE_DEVICE_AND_HOST static residue_type next_value_mod_2(residue_type v_n_mod_p, index_type n);
        /// Computes the shift by which Sm3(n) needs to be multiplied in order to append n + 1
        /// e.g. Sm3(8) * shift + 9 = Sm3(9), shift is 3, so the returned value is simply
        /// ceil(log3(n + 1)).
        CUTRIALDIVE_DEVICE_AND_HOST static index_type compute_shift(index_type n);
    };
}

namespace cutrialdive {

    // @todo: commonalities with base 10 below
    // especially next_value_mod and next_value_mod_mu
    // 

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    typename number_sequence<num_seq_id::smarandache, 3>::index_type
    number_sequence<num_seq_id::smarandache, 3>::compute_shift(index_type n)
    {
        auto m = n + 1;
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
        index_type result = 3;
        while(m /= 3) {
            result *= 3;
        }
        return result;
    }


    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    typename number_sequence<num_seq_id::smarandache, 3>::residue_type
    number_sequence<num_seq_id::smarandache, 3>::value_mod_2(index_type n)
    {
        // Sm3(n) mod 2 = 1 if and only if n mod 4 belongs to {1, 2}
        return (n & 3) && ((n & 3) != 3);
    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    typename number_sequence<num_seq_id::smarandache, 3>::residue_type
    number_sequence<num_seq_id::smarandache, 3>::next_value_mod(
        residue_type v_n_mod_p,
        index_type n,
        residue_type p
    )
    {
        // v_n_mod_p is Sm3(n - 1) mod p, need to return Sm3(n) mod p
        // Sm3(n + 1) = Sm3(n) * 3^(ceil(log3(n + 1))) + n + 1
        return (__uint128_t(v_n_mod_p) * compute_shift(n) + n + 1) % p;
    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    typename number_sequence<num_seq_id::smarandache, 3>::residue_type
    number_sequence<num_seq_id::smarandache, 3>::next_value_mod_mu(
        residue_type v_n_mod_p,
        index_type n,
        residue_type p,
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
    typename number_sequence<num_seq_id::smarandache, 3>::residue_type
    number_sequence<num_seq_id::smarandache, 3>::next_value_mod_2(residue_type v_n_mod_p, index_type n)
    {
        // v_n_mod_p is guaranteed to be a residue mod 2
        // the shift is a power of 3 thus always odd
        return v_n_mod_p ^ (n & 1) ^ 1;
    }
}
