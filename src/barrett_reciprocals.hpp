/*
* Created on 2026.05.13
* Copyright (c) Youcef Lemsafer
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
        if constexpr(std::is_same_v<BarrettMuType, barrett_mu64_t>) {
            return barrett_mu64_t{mu64(divisor)};
        } else if constexpr(std::is_same_v<BarrettMuType, barrett_mu128_t>) {
            return barrett_mu128_t{mu128(divisor)};
        } else if constexpr(std::is_same_v<BarrettMuType, barrett_mu_both_t>) {
            return barrett_mu_both_t{mu128(divisor), mu64(divisor)};
        }
    }

    template <typename NumberSequenceT>
    CUTRIALDIVE_DEVICE CUTRIALDIVE_INLINE
    auto compute_barrett_mu_numseq(uint64_t divisor)
    {
        if constexpr(HasValueModMuImpl<NumberSequenceT, barrett_mu64_t>
                    || HasNextValueModMuImpl<NumberSequenceT, barrett_mu64_t>) {
            return compute_barrett_mu<barrett_mu64_t>(divisor);
        } else if constexpr(HasValueModMuImpl<NumberSequenceT, barrett_mu128_t>
                           || HasNextValueModMuImpl<NumberSequenceT, barrett_mu128_t>) {
            return compute_barrett_mu<barrett_mu128_t>(divisor);
        } else if constexpr(HasValueModMuImpl<NumberSequenceT, barrett_mu_both_t>
                            || HasNextValueModMuImpl<NumberSequenceT, barrett_mu_both_t>) {
            return compute_barrett_mu<barrett_mu_both_t>(divisor);
        } else {
            return no_barrett_t{};
        }
    }
}
