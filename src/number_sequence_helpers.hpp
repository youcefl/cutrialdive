/*
* MIT License
* Created on 2026.06.03
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#pragma once

#include "common_defines.h"
#include "number_sequence.hpp"

namespace cutrialdive {

    template <typename NumberSequenceT>
    CUTRIALDIVE_DEVICE_AND_HOST CUTRIALDIVE_INLINE
    auto value_mod_2(NumberSequenceT seq,  typename NumberSequenceT::index_type n)
    {
        if constexpr(HasValueModTwo<NumberSequenceT>) {
            return seq.value_mod_2(n);
        } else {
            return seq.value_mod(n, 2);
        }
    }

    template <typename NumberSequenceT>
    CUTRIALDIVE_DEVICE_AND_HOST CUTRIALDIVE_INLINE
    auto next_value_mod_2(NumberSequenceT seq, typename NumberSequenceT::residue_type r, typename NumberSequenceT::index_type n)
    {
        if constexpr(HasNextValueModTwo<NumberSequenceT>) {
            return seq.next_value_mod_2(r, n);
        } else {
            return seq.next_value_mod(r, n, 2);
        }
    }

    template <typename NumberSequenceT>
    auto short_name(NumberSequenceT seq)
    {
        if constexpr(HasShortName<NumberSequenceT>) {
            return seq.short_name();
        } else {
            return get_math_sequence(seq).short_name();
        }
    }
}

