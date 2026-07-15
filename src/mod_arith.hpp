/*
* MIT License
* Created on 2026.07.15
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#pragma once

#include "common_defines.h"
#include "modular_arithmetic_detail.hpp"


namespace cutrialdive {

    struct mod_arith
    {
        using value_type = uint64_t;
        uint64_t d;
        mu_both_t mu;
        CUTRIALDIVE_DEVICE_AND_HOST value_type mul_by_2(value_type a) const;
        CUTRIALDIVE_DEVICE_AND_HOST value_type add(value_type a, value_type b) const;
        CUTRIALDIVE_DEVICE_AND_HOST value_type sub(value_type a, value_type b) const;
        CUTRIALDIVE_DEVICE_AND_HOST value_type mul(value_type a, value_type b) const;
        CUTRIALDIVE_DEVICE_AND_HOST value_type sqr(value_type a) const;
    };

}

namespace cutrialdive {

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    mod_arith::value_type mod_arith::mul_by_2(value_type a) const
    {
        return doublemod(a, d);
    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    mod_arith::value_type mod_arith::add(value_type a, value_type b) const
    {
        return addmod(a, b, d);
    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    mod_arith::value_type mod_arith::sub(value_type a, value_type b) const
    {
        return submod(a, b, d);
    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    mod_arith::value_type mod_arith::mul(value_type a, value_type b) const
    {
        return mulmod(a, b, d, mu.mu64, mu.mu128);
    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    mod_arith::value_type mod_arith::sqr(value_type a) const
    {
        return mulmod(a, a, d, mu.mu64, mu.mu128);
    }

}
