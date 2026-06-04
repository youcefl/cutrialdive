/*
* Creation date: 2026.06.02
* Created by Youcef Lemsafer
*/
#pragma once

#include <concepts>

#include "barrett_mu_types.hpp"


namespace cutrialdive {

    /// N.B.: below S(n) is used to designate the sequence being modelled by Seq.

    /// Boolean indicating whether S(n0) mod d has to be computed by first computing S(n0) = seq.value(n0)
    /// and then reducing that value modulo d. If false then S(n0) mod d is computed by calling
    /// seq.value_mod(n0, d) or seq.value_mod_mu(n0, d).
    template <typename Seq>
    concept InitializeFromValue = Seq::initialize_from_value;

    /// True if and only if the sequence has a function n -> S(n)
    template <typename Seq>
    concept HasValue =
        requires(
            Seq seq,
            typename Seq::index_type n
        )
    {
        {
            seq.value(n)
        } -> std::same_as<typename Seq::value_type>;
    };

    /// True if and only if the sequence has a function (n, d) -> S(n) mod d
    /// If true the function will be used to compute S(n0) mod d unless
    /// HasValueModMu is true.
    template <typename Seq>
    concept HasValueMod =
        requires(
            Seq seq,
            typename Seq::index_type n,
            typename Seq::residue_type d
        )
    {
        {
            seq.value_mod(n, d)
        } -> std::same_as<typename Seq::residue_type>;
    };

    /// True if and only if the sequence has a function (r, n, d) -> S(n + 1) mod d
    /// where r is S(n) mod d.
    template <typename Seq>
    concept HasNextValueMod =
        requires(
            Seq seq,
            typename Seq::residue_type r,
            typename Seq::index_type n,
            typename Seq::residue_type d
        )
    {
        {
            seq.next_value_mod(r, n, d)
        } -> std::same_as<typename Seq::residue_type>;
    };

    /// True if and only if the sequence has a function n -> S(n) mod 2
    /// If true then the engine uses the function for evaluating S(n) mod 2
    template <typename Seq>
    concept HasValueModTwo =
        requires(
            Seq seq,
            typename Seq::index_type n
        )
    {
        {
            // Returns S(n) mod 2, given n and d
            seq.value_mod_2(n)
        } -> std::same_as<typename Seq::residue_type>;
    };

    /// True if and only if the sequence has a function (r, n) -> S(n + 1) mod 2 where r = S(n) mod 2
    /// If true then the engine uses the function for propagating residue mod 2
    template <typename Seq>
    concept HasNextValueModTwo =
        requires(
            Seq seq,
            typename Seq::residue_type r,
            typename Seq::index_type n
        )
    {
        {
            // Returns S(n+1) mod 2, given r = S(n) mod 2 and n
            seq.next_value_mod_2(r, n)
        } -> std::same_as<typename Seq::residue_type>;
    };

    template <typename Seq, typename BarrettMuT>
    concept HasValueModMuImpl =
        requires(
            Seq seq,
            typename Seq::index_type n,
            typename Seq::residue_type d,
            BarrettMuT mu
        )
    {
        {
            seq.value_mod_mu(n, d, mu)
        } -> std::same_as<typename Seq::residue_type>;
    };

    /// True if and only if the sequence has a function (n, d, mu) -> S(n) mod d where mu is the Barrett reciprocal of d
    /// If true the engine will use this function to evaluate S(n0) mod d
    template <typename Seq>
    concept HasValueModMu =
           HasValueModMuImpl<Seq, mu64_t>
        || HasValueModMuImpl<Seq, mu128_t>
        || HasValueModMuImpl<Seq, mu_both_t>
    ;

    template <typename Seq, typename BarrettMuT>
    concept HasNextValueModMuImpl =
        requires(
            Seq seq,
            typename Seq::residue_type r,
            typename Seq::index_type n,
            typename Seq::residue_type d,
            BarrettMuT mu
        )
    {
        {
            seq.next_value_mod_mu(r, n, d, mu)
        } -> std::same_as<typename Seq::residue_type>;
    };

    /// True if and only if the sequence has a function (r, n, d, mu) -> S(n) mod d where r = S(n) mod d and
    ///  mu is the Barrett reciprocal of d.
    /// If true the engine will use this function to evaluate S(n + 1) mod d given S(n) mod d
    template <typename Seq>
    concept HasNextValueModMu =
           HasNextValueModMuImpl<Seq, mu64_t>
        || HasNextValueModMuImpl<Seq, mu128_t>
        || HasNextValueModMuImpl<Seq, mu_both_t>
    ;


    template <typename Seq>
    concept HaveConsistentMuType =
        !HasValueModMu<Seq> || !HasNextValueModMu<Seq>
        || ((HasValueModMuImpl<Seq, mu64_t> == HasNextValueModMuImpl<Seq, mu64_t>)
            && (HasValueModMuImpl<Seq, mu128_t> == HasNextValueModMuImpl<Seq, mu128_t>)
            && (HasValueModMuImpl<Seq, mu_both_t> == HasNextValueModMuImpl<Seq, mu_both_t>)
        )
    ;

    /// True if and only if the sequence has a function print_expression(n, o) where n is the index
    /// and o an output stream. The function writes an expression which evaluates to S(n) into o and returns o
    template <typename Seq>
    concept HasExpressionPrinter =
        requires(
            Seq seq,
            typename Seq::index_type n,
            std::ostream & out
        )
    {
        {
            seq.print_expression(n, out)
        } -> std::same_as<std::ostream &>;
    };

    /// True if and only if the sequence has a function print_value(n, o) where n is the index
    /// and o an output stream. The function writes S(n) in base 10 into o and returns o
    template <typename Seq>
    concept HasValuePrinter =
        requires(
            Seq seq,
            typename Seq::index_type n,
            std::ostream & out
        )
    {
        {
            seq.print_value(n, out)
        } -> std::same_as<std::ostream &>;
    };

}
