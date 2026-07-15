/*
* MIT License
* Created on 2026.06.02
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#pragma once

#include <concepts>
#include <ostream>
#include <string_view>

#include "barrett_mu_types.hpp"


namespace cutrialdive {

    /// N.B.: below S(n) is used to designate the sequence being modelled by Seq.

    template <typename T>
    inline
    decltype(auto) get_math_sequence(T&& seq)
    {
        if constexpr(requires { seq.math_sequence(); }) {
            return seq.math_sequence();
        } else {
            return std::forward<T>(seq);
        }
    }

    namespace details {

        template <typename Seq>
        using math_sequence_type = typename std::remove_cvref_t<
                decltype(get_math_sequence(std::declval<Seq>()))
            >;

    }

    /// Struct used when sequence has no state
    struct no_state_t {};

    /// True if and only if the sequence is stateful
    template <typename Seq>
    concept HasState = Seq::has_state;

    template <typename Seq>
    concept HasMakeState =
        requires(
            Seq seq,
            typename Seq::index_type n,
            typename Seq::residue_type d
        )
    {
        {
            seq.make_state(n, d)
        } -> std::same_as<typename Seq::state_type>;
    };
    template <typename Seq>
    concept HasMakeStateWithMu =
        requires(
            Seq seq,
            typename Seq::index_type n,
            typename Seq::residue_type d,
            typename Seq::mu_type mu
        )
    {
        {
            seq.make_state(n, d, mu)
        } -> std::same_as<typename Seq::state_type>;
    };



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
    template <typename Seq>
    concept HasValueModWithState =
        requires(
            Seq seq,
            typename Seq::index_type n,
            typename Seq::residue_type d,
            typename Seq::state_type & state
        )
    {
        {
            seq.value_mod(n, d, state)
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
    template <typename Seq>
    concept HasNextValueModWithState =
        requires(
            Seq seq,
            typename Seq::residue_type r,
            typename Seq::index_type n,
            typename Seq::residue_type d,
            typename Seq::state_type & state
        )
    {
        {
            seq.next_value_mod(r, n, d, state)
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
            // Returns S(n) mod 2, given n 
            seq.value_mod_2(n)
        } -> std::same_as<typename Seq::residue_type>;
    };
    template <typename Seq>
    concept HasValueModTwoWithState =
        requires(
            Seq seq,
            typename Seq::index_type n,
            typename Seq::state_type & state
        )
    {
        {
            // Returns S(n) mod 2, given n
            seq.value_mod_2(n, state)
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
    template <typename Seq>
    concept HasNextValueModTwoWithState =
        requires(
            Seq seq,
            typename Seq::residue_type r,
            typename Seq::index_type n,
            typename Seq::state_type & state
        )
    {
        {
            // Returns S(n+1) mod 2, given r = S(n) mod 2 and n
            seq.next_value_mod_2(r, n, state)
        } -> std::same_as<typename Seq::residue_type>;
    };

    /// True if and only if the sequence has a function (n, d, mu) -> S(n) mod d where mu is the Barrett reciprocal of d
    /// If true the engine will use this function to evaluate S(n0) mod d
    template <typename Seq>
    concept HasValueModMu =
        requires(
            Seq seq,
            typename Seq::index_type n,
            typename Seq::residue_type d,
            typename Seq::mu_type mu
        )
    {
        {
            seq.value_mod_mu(n, d, mu)
        } -> std::same_as<typename Seq::residue_type>;
    };
    /// Same as HasValueModMu for stateful sequences
    template <typename Seq>
    concept HasValueModMuWithState =
        requires(
            Seq seq,
            typename Seq::index_type n,
            typename Seq::residue_type d,
            typename Seq::mu_type mu,
            typename Seq::state_type & state
        )
    {
        {
            seq.value_mod_mu(n, d, mu, state)
        } -> std::same_as<typename Seq::residue_type>;
    };

    /// True if and only if the sequence has a function (r, n, d, mu) -> S(n) mod d where r = S(n) mod d and
    ///  mu is the Barrett reciprocal of d.
    /// If true the engine will use this function to evaluate S(n + 1) mod d given S(n) mod d
    template <typename Seq>
    concept HasNextValueModMu =
        requires(
            Seq seq,
            typename Seq::residue_type r,
            typename Seq::index_type n,
            typename Seq::residue_type d,
            typename Seq::mu_type mu
        )
    {
        {
            seq.next_value_mod_mu(r, n, d, mu)
        } -> std::same_as<typename Seq::residue_type>;
    };
    /// Same as HasNextValueModMu for stateful sequences
    template <typename Seq>
    concept HasNextValueModMuWithState =
        requires(
            Seq seq,
            typename Seq::residue_type r,
            typename Seq::index_type n,
            typename Seq::residue_type d,
            typename Seq::mu_type mu,
            typename Seq::state_type & state
        )
    {
        {
            seq.next_value_mod_mu(r, n, d, mu, state)
        } -> std::same_as<typename Seq::residue_type>;
    };


    /// True if and only if the sequence has a function print_expression(n, o) where n is the index
    /// and o an output stream. The function writes an expression which evaluates to S(n) into o and returns o
    template <typename Seq>
    concept HasExpressionPrinter =
        requires(
            Seq seq,
            typename details::math_sequence_type<Seq>::index_type n,
            std::ostream & out
        )
    {
        {
            get_math_sequence(seq).print_expression(n, out)
        } -> std::same_as<std::ostream &>;
    };

    /// True if and only if the sequence has a function print_value(n, o) where n is the index
    /// and o an output stream. The function writes S(n) in base 10 into o and returns o
    template <typename Seq>
    concept HasValuePrinter =
        requires(
            Seq seq,
            typename details::math_sequence_type<Seq>::index_type n,
            std::ostream & out
        )
    {
        {
            get_math_sequence(seq).print_value(n, out)
        } -> std::same_as<std::ostream &>;
    };

    template <typename Seq>
    concept HasShortName =
        requires(
            Seq seq
        )
    {
        {
            seq.short_name()
        } -> std::convertible_to<std::string_view>;
    };

    template <typename Seq>
    concept HasSequenceShortName =
        HasShortName<Seq>
     || HasShortName<details::math_sequence_type<Seq>>;

}
