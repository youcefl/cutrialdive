/*
* Creation date: 2026.04.29
* Created by Youcef Lemsafer
*/
#pragma once

#include <ostream>
#include <string_view>
#include "number_sequence_detail.hpp"

namespace cutrialdive {

    template <typename Seq>
    concept NumberSequence = requires(Seq seq)
    {
        typename Seq::index_type;
        typename Seq::residue_type;

        // A short name for the sequence e.g. M for Mersenne, Sm for Smarandache.
        { seq.short_name() }
            -> std::convertible_to<std::string_view>;
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

}


