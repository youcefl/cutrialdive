/*
* Creation date: 2026.04.29
* Created by Youcef Lemsafer
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
    )
    &&
    (
        // Need to be able to compute residue mod 2
        (InitializeFromValue<Seq> && HasValue<Seq>)
        || HasValueMod<Seq>
        || HasValueModTwo<Seq>
    )
    &&
    (
        // Needed for propagation of residues
        HasNextValueMod<Seq> || HasNextValueModMu<Seq>
    )
    &&
    (
        // Need to be able to propagate residue mod 2
        HasNextValueMod<Seq> || HasNextValueModTwo<Seq>
    )
        // If both value_mod_mu and next_value_mod_mu are provided they have to have the same mu type
    &&  HaveConsistentMuType<Seq>
    ;

    template <typename Seq>
    concept NumberSequence = 
        PureMathSequence<details::math_sequence_type<Seq>>
     && HasSequenceShortName<Seq>;
}

