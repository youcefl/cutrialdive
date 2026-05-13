/*
* Created on 2026.05.13
* Copyright (c) Youcef Lemsafer
*/
#pragma once

#include "common_defines.h"
#include "barrett_mu_types.hpp"

namespace cutrialdive {

    /// Returns the Barrett reciprocal of @param divisor as required by the number sequence
    /// @pre divisor is odd
    template <typename NumberSequenceT>
    CUTRIALDIVE_DEVICE CUTRIALDIVE_INLINE
    auto compute_barrett_mu(uint64_t divisor)
    {
        using MuType = typename NumberSequenceT::barrett_mu_type;
        if constexpr(std::is_same_v<MuType, barrett_mu64_t>) {
            return barrett_mu64_t{mu64(divisor)};
        } else if constexpr(std::is_same_v<MuType, barrett_mu128_t>) {
            return barrett_mu128_t{mu128(divisor)};
        } else if constexpr(std::is_same_v<MuType, barrett_mu_both_t>) {
            return barrett_mu_both_t{mu128(divisor), mu64(divisor)};
        } else if constexpr(std::is_same_v<MuType, no_barrett_t>) {
            return no_barrett_t{};
        }
    }
}
