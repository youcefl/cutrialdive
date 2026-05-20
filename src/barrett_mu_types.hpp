/*
* Creation date: 2026.05.13
* Created by Youcef Lemsafer
*/
#pragma once

#include <cstdint>

namespace cutrialdive {

    /// No Barrett multiplication
    struct no_barrett_t {};

    /// The Barrett 64-bit reciprocal
    using barrett_mu64_t = uint64_t;

    /// The Barrett 128-bit reciprocal
    using barrett_mu128_t = __uint128_t;

    /// Both the 64-bit and 128-bit reciprocals
    struct barrett_mu_both_t
    {
        barrett_mu128_t mu128;
        uint64_t mu64;
    };

}
