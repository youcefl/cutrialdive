/*
* Creation date: 2026.04.29
* Created by Youcef Lemsafer
*/
#pragma once

#include "mode.hpp"

namespace cutrialdive {

    template <typename Seq>
    concept NumberSequence = requires {
        typename Seq::index_type;
        typename Seq::value_type;
        typename Seq::residue_type;
        { Seq::value(std::declval<typename Seq::index_type>()) } -> std::same_as<typename Seq::value_type>;
        { Seq::next_value_mod_mu(std::declval<typename Seq::residue_type>(),
            std::declval<typename Seq::index_type>(),
            std::declval<typename Seq::residue_type>(),
        ) } -> std::same_as<typename Seq::residue_type>;
    };

    template <mode_flag flag>
    struct number_sequence;
}

