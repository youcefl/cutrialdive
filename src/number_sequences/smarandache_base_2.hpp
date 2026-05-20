/*
 * Creation date: 2026.05.15
 * Created by Youcef Lemsafer
 */
#pragma once

#include <cstdint>
#include <ostream>
#if !CUTRIALDIVE_IS_CUDA
#include <bit>
#endif
#include "common_defines.h"
#include "barrett_mu_types.hpp"
#include "number_sequence.hpp"
#include "hgint.hpp"
#include "modular_arithmetic_detail.hpp"

namespace cutrialdive {

    template <>
    struct number_sequence<num_seq_id::smarandache, 2>
    {
        static constexpr uint32_t base = 2;
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
    };
}

namespace cutrialdive {

#if CUTRIALDIVE_IS_CUDA
#  define CUTRIALDIVE_BIT_WIDTH __clzll
#else
#  define CUTRIALDIVE_BIT_WIDTH std::bit_width
#endif

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    typename number_sequence<num_seq_id::smarandache, 2>::residue_type
    number_sequence<num_seq_id::smarandache, 2>::value_mod_2(index_type n)
    {
        return n & 1;
    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    typename number_sequence<num_seq_id::smarandache, 2>::residue_type
    number_sequence<num_seq_id::smarandache, 2>::next_value_mod(residue_type v_n_mod_p, index_type n, residue_type p)
    {
        auto shift = sizeof(n) * 8 - CUTRIALDIVE_BIT_WIDTH(n + 1);
        return ((__uint128_t(v_n_mod_p) << shift) + n + 1) % p;
    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    typename number_sequence<num_seq_id::smarandache, 2>::residue_type
    number_sequence<num_seq_id::smarandache, 2>::next_value_mod_mu(residue_type v_n_mod_p, index_type n, residue_type p, barrett_mu_type mu_p)
    {
        auto shift = sizeof(n) * 8 - CUTRIALDIVE_BIT_WIDTH(n + 1);
        auto a = (__uint128_t(v_n_mod_p) << shift) + n + 1;
        if(a >> 64) {
            auto q = mul128hi(a, mu_p.mu128);
            auto r = a - q * p;
            return r >= p ? r - p : r;
        }
        auto q = (__uint128_t(a) * mu_p.mu64) >> 64;
        auto r = a - q * p;
        return r >= p ? r - p : r;
    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    typename number_sequence<num_seq_id::smarandache, 2>::residue_type
    number_sequence<num_seq_id::smarandache, 2>::next_value_mod_2(residue_type v_n_mod_p, index_type n)
    {
        return (n & 1) ^ 1;
    }
}
