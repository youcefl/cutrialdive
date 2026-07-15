/*
* MIT License
* Creation date: 2026.05.13
* Created by Youcef Lemsafer
* See LICENSE file for details.
*/
#pragma once

#include <cstdint>
#include "common_defines.h"

namespace cutrialdive {

    /// No Barrett multiplication
    struct no_mu_t {};

    /// The Barrett 64-bit reciprocal
    struct mu64_t {
        uint64_t v;
        CUTRIALDIVE_DEVICE_AND_HOST
        explicit mu64_t(uint64_t x) : v(x) {}
    };

    /// The Barrett 128-bit reciprocal
    struct mu128_t {
        __uint128_t v;
        CUTRIALDIVE_DEVICE_AND_HOST
        explicit mu128_t(__uint128_t x) : v(x) {}
    };

    /// Both the 64-bit and 128-bit reciprocals
    struct mu_both_t
    {
        __uint128_t mu128;
        uint64_t mu64;
    };

}
