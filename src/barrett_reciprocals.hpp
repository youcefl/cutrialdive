/*
* MIT License
* Created on 2026.05.13
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#pragma once

#include "common_defines.h"
#include "barrett_mu_types.hpp"
#include "modular_arithmetic_detail.hpp"
#include "number_sequence.hpp"

namespace cutrialdive {

    /// Returns the Barrett reciprocal of @param divisor
    /// @pre divisor is odd
    template <typename BarrettMuType>
    CUTRIALDIVE_DEVICE CUTRIALDIVE_INLINE
    auto compute_barrett_mu(uint64_t divisor)
    {
        if constexpr(std::is_same_v<BarrettMuType, mu64_t>) {
            return mu64_t{mu64(divisor)};
        } else if constexpr(std::is_same_v<BarrettMuType, mu128_t>) {
            return mu128_t{mu128(divisor)};
        } else if constexpr(std::is_same_v<BarrettMuType, mu_both_t>) {
            return mu_both_t{mu128(divisor), mu64(divisor)};
        }
    }

    template <typename NumberSequenceT>
    CUTRIALDIVE_DEVICE CUTRIALDIVE_INLINE
    auto compute_barrett_mu_numseq(uint64_t divisor)
    {
        if constexpr (requires { typename NumberSequenceT::mu_type; }) {
            return compute_barrett_mu<typename NumberSequenceT::mu_type>(divisor);
        } else {
            return no_mu_t{};
        }
    }
}
