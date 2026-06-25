/*
* MIT License
* Created on 2026.06.11
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/

#include <catch2/catch_test_macros.hpp>

#include "number_sequences/proth.hpp"
#include "barrett_reciprocals.hpp"


namespace ctd = cutrialdive;

TEST_CASE("Basic tests for Proth numbers 1")
{
    auto k123 = ctd::proth{123};

    REQUIRE(k123.short_name() == std::string{"Pr_123"});

    auto k123Seq = k123.math_sequence();

    REQUIRE(k123Seq.value(0) == 124);
    REQUIRE(k123Seq.value(1) == 247);
    REQUIRE(k123Seq.value(2) == 493);
}

TEST_CASE("Basic tests for Proth numbers 2")
{
    auto k78557Seq = ctd::proth{78557}.math_sequence();

    REQUIRE(k78557Seq.value(0) == 78558);
    REQUIRE(k78557Seq.value_mod_2(0) == 0);
    REQUIRE(k78557Seq.value_mod_2(1) == 1);
    REQUIRE(k78557Seq.value_mod_2(2) == 1);
    REQUIRE(k78557Seq.value_mod_2(1276481311) == 1);

    auto bm3 = ctd::compute_barrett_mu<ctd::mu_both_t>(3);
    REQUIRE(k78557Seq.value_mod_mu(0, 3, bm3) == 78558 % 3);
    REQUIRE(k78557Seq.value_mod_mu(1, 3, bm3) == 2);
    REQUIRE(k78557Seq.value_mod_mu(2, 3, bm3) == 0);
    REQUIRE(k78557Seq.value_mod_mu(3, 3, bm3) == 2);
    REQUIRE(k78557Seq.value_mod_mu(797827, 3, bm3) == 2);
    REQUIRE(k78557Seq.value_mod_mu(877532, 3, bm3) == 0);
}

TEST_CASE("Basic tests for Proth numbers 3")
{
    auto k271129Seq = ctd::proth{271129}.math_sequence();

    auto bm5 = ctd::compute_barrett_mu<ctd::mu_both_t>(5);
    REQUIRE(k271129Seq.value_mod_mu(0, 5, bm5) == 271130 % 5);
    REQUIRE(k271129Seq.value_mod_mu(1, 5, bm5) == 4);
    REQUIRE(k271129Seq.value_mod_mu(2, 5, bm5) == 2);
    REQUIRE(k271129Seq.value_mod_mu(3, 5, bm5) == 3);
    REQUIRE(k271129Seq.value_mod_mu(4, 5, bm5) == 0);
    REQUIRE(k271129Seq.value_mod_mu(6700224, 5, bm5) == 0);
    REQUIRE(k271129Seq.value_mod_mu(4924567891ull, 5, bm5) == 3);
}

TEST_CASE("Modular recurrence tests for Proth numbers 1")
{
    auto k271129Seq = ctd::proth{271129}.math_sequence();

    auto bm5 = ctd::compute_barrett_mu<ctd::mu_both_t>(5);
    auto n0 = 6700224u;
    auto r0 = k271129Seq.value_mod_mu(n0, 5, bm5);
    REQUIRE(k271129Seq.next_value_mod(r0, n0, 5) == 4);
    REQUIRE(k271129Seq.next_value_mod(4, n0 + 1, 5) == 2);
    REQUIRE(k271129Seq.next_value_mod(2, n0 + 2, 5) == 3);
    REQUIRE(k271129Seq.next_value_mod(3, n0 + 3, 5) == 0);
    REQUIRE(k271129Seq.next_value_mod(0, n0 + 4, 5) == 4);
    REQUIRE(k271129Seq.next_value_mod(4, n0 + 5, 5) == 2);
    REQUIRE(k271129Seq.next_value_mod(2, n0 + 6, 5) == 3);
    REQUIRE(k271129Seq.next_value_mod(3, n0 + 7, 5) == 0);
}
