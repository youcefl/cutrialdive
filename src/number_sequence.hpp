/*
* Creation date: 2026.04.29
* Created by Youcef Lemsafer
*/
#pragma once

#include <ostream>
#include "num_seq_id.hpp"

namespace cutrialdive {

    template <typename Seq>
    concept NumberSequence = requires {
        typename Seq::index_type;
        typename Seq::value_type;
        typename Seq::residue_type;
        typename Seq::barrett_mu_type;
        
        // A short name for the sequence e.g. M for Mersenne, Sm for Smarandache.
        { Seq::short_name() } 
            -> std::same_as<char const*>;

        // Prints the value of S(n) in base 10
        { Seq::print_value(std::declval<typename Seq::index_type>(), std::declval<std::ostream&>()) } 
            -> std::same_as<void>;

        // Prints an expression that evaluates to S(n)
        { Seq::print_expression(std::declval<typename Seq::index_type>(), std::declval<std::ostream&>()) } 
            -> std::same_as<void>;
        
        // Returns S(n)
        { Seq::value(std::declval<typename Seq::index_type>()) }
            -> std::same_as<typename Seq::value_type>;

        // Returns S(n) mod 2
        { Seq::value_mod_2(std::declval<typename Seq::index_type>()) }
            -> std::same_as<typename Seq::residue_type>;

        // Given S(n) mod 2 return S(n+1) mod 2
        { Seq::next_value_mod_2(
            std::declval<typename Seq::residue_type>(),
            std::declval<typename Seq::index_type>()
        ) } -> std::same_as<typename Seq::residue_type>;
        
        // Returns S(n+1) mod d given S(n) mod d, n and d
        // d is guaranteed to be odd
        { Seq::next_value_mod(
            std::declval<typename Seq::residue_type>(), // S(n) mod d
            std::declval<typename Seq::index_type>(),   // n
            std::declval<typename Seq::residue_type>()  // d
        ) } -> std::same_as<typename Seq::residue_type>;
        
        // Modular recurrence with Barrett multiplier
        // Given S(n) mod d, n, d and mu(d) return S(n+1) mod d
        // d is guaranteed to be odd
        // mu(d) is the Barrett reciprocal of d
        { Seq::next_value_mod_mu(
            std::declval<typename Seq::residue_type>(), // S(n) mod d
            std::declval<typename Seq::index_type>(),   // n
            std::declval<typename Seq::residue_type>(), // d
            std::declval<typename Seq::barrett_mu_type>()  // mu(d)
        ) } -> std::same_as<typename Seq::residue_type>;
    };
}

