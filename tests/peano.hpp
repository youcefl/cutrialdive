/*
* Created on 2026.06.04
* Copyright (c) Youcef Lemsafer
*/
#pragma once

#include <cstdint>
#include "common_defines.h"
#include "modular_arithmetic_detail.hpp"
#include "barrett_mu_types.hpp"

namespace cutrialdive::tests {

    enum class function_to_define
    {
        value,
        value_mod,
        value_mod_2,
        value_mod_mu64,
        value_mod_mu128,
        next_value_mod,
        next_value_mod_2,
        next_value_mod_mu128
    };

    template <function_to_define func>
    struct peano_function {};

    template <>
    struct peano_function<function_to_define::value>
    {
        static uint64_t value(uint64_t n)
        {
            return n;
        }
    };

    template <>
    struct peano_function<function_to_define::value_mod>
    {
        CUTRIALDIVE_DEVICE_AND_HOST
        static uint64_t value_mod(uint64_t n, uint64_t d)
        {
            return n % d;
        }
    };

    template <>
    struct peano_function<function_to_define::value_mod_2>
    {
        CUTRIALDIVE_DEVICE_AND_HOST
        static uint64_t value_mod_2(uint64_t n)
        {
            return n & 1;
        }
    };

    template <>
    struct peano_function<function_to_define::value_mod_mu64>
    {
        CUTRIALDIVE_DEVICE_AND_HOST
        static uint64_t value_mod_mu(uint64_t n, uint64_t d, mu64_t mu)
        {
            auto q = (__uint128_t(n) * mu.v) >> 64;
            auto r = n - q * d;
            return r >= d ? r - d : r;
        }
    };

    template <>
    struct peano_function<function_to_define::value_mod_mu128>
    {
        CUTRIALDIVE_DEVICE_AND_HOST
        static uint64_t value_mod_mu(uint64_t n, uint64_t d, mu128_t mu)
        {
            auto q = mul128hi(n, mu.v);
            auto r = __uint128_t(n) - q * __uint128_t(d);
            return uint64_t(r >= d ? r - d : r);
        }
    };

    template <>
    struct peano_function<function_to_define::next_value_mod>
    {
        CUTRIALDIVE_DEVICE_AND_HOST
        static uint64_t next_value_mod(uint64_t r, uint64_t, uint64_t d)
        {
            r = r + 1;
            return r == d ? 0 : r; 
        }
    };

    template <>
    struct peano_function<function_to_define::next_value_mod_2>
    {
        CUTRIALDIVE_DEVICE_AND_HOST
        static uint64_t next_value_mod_2(uint64_t r, uint64_t n)
        {
            return r == 0 ? 1 : 0;
        }
    };

    template <>
    struct peano_function<function_to_define::next_value_mod_mu128>
    {
        CUTRIALDIVE_DEVICE_AND_HOST
        static uint64_t next_value_mod_mu(uint64_t r, uint64_t n, uint64_t d, mu128_t mu)
        {
            r = r + 1;
            auto q = mul128hi(r, mu.v);
            auto s = r - q * d;
            return uint64_t(s >= d ? s - d : s);
        }
    };

    template <function_to_define... funcs>
    struct peano : peano_function<funcs>...
    {
        using index_type = uint64_t;
        using residue_type = uint64_t;

        static char const* short_name()
        {
            return "P";
        }
    };

}
