/*
* MIT License
* Created on 2026.04.29
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#pragma once

#include <ostream>
#include <string_view>
#include <type_traits>
#include <utility>

#include "number_sequence_detail.hpp"


namespace cutrialdive {

    template <typename Seq>
    concept PureMathSequence = requires(Seq seq)
    {
        typename Seq::index_type;
        typename Seq::residue_type;
    }
    &&
    (
        // Need one of these for initialization
        (InitializeFromValue<Seq> && HasValue<Seq>)
        || HasValueMod<Seq>
        || HasValueModMu<Seq>
        || (HasState<Seq> && (HasValueModWithState<Seq> || HasValueModMuWithState<Seq>))
    )
    &&
    (
        // Need to be able to compute residue mod 2
        (InitializeFromValue<Seq> && HasValue<Seq>)
        || HasValueMod<Seq>
        || HasValueModTwo<Seq>
        || (HasState<Seq> && (HasValueModWithState<Seq> || HasValueModTwoWithState<Seq>))
    )
    &&
    (
        // Needed for propagation of residues
        HasNextValueMod<Seq> || HasNextValueModMu<Seq>
        || (HasState<Seq> && (HasNextValueModWithState<Seq> || HasNextValueModMuWithState<Seq>))
    )
    &&
    (
        // Need to be able to propagate residue mod 2
        HasNextValueMod<Seq> || HasNextValueModTwo<Seq>
        || (HasState<Seq> && (HasNextValueModWithState<Seq> || HasNextValueModTwoWithState<Seq>))
    )
    ;

    template <typename Seq>
    concept NumberSequence = 
        PureMathSequence<details::math_sequence_type<Seq>>
     && HasSequenceShortName<Seq>;
}

