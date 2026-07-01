/*
* MIT License
* Created on 2026.06.28
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/

#include <catch2/catch_test_macros.hpp>

#include "number_sequences/riesel_proth.hpp"
#include "barrett_reciprocals.hpp"


namespace ctd = cutrialdive;

TEST_CASE("Basic tests for Riesel sequence 1")
{
    auto k37 = ctd::riesel{37};

    REQUIRE(k37.short_name() == std::string{"R_37"});

    auto k37Seq = k37.math_sequence();

    REQUIRE(k37Seq.value(0) == 36);
    REQUIRE(k37Seq.value(1) == 73);
    REQUIRE(k37Seq.value(2) == 147);
}

TEST_CASE("Basic tests for Riesel sequence 2")
{
    auto k78557Seq = ctd::riesel{78557}.math_sequence();

    REQUIRE(k78557Seq.value(0) == 78556);
    REQUIRE(k78557Seq.value_mod_2(0) == 0);
    REQUIRE(k78557Seq.value_mod_2(1) == 1);
    REQUIRE(k78557Seq.value_mod_2(2) == 1);
    REQUIRE(k78557Seq.value_mod_2(1276481311) == 1);

    auto bm3 = ctd::compute_barrett_mu<ctd::mu_both_t>(3);
    REQUIRE(k78557Seq.value_mod_mu(0, 3, bm3) == 78556 % 3);
    REQUIRE(k78557Seq.value_mod_mu(1, 3, bm3) == 0);
    REQUIRE(k78557Seq.value_mod_mu(2, 3, bm3) == 1);
    REQUIRE(k78557Seq.value_mod_mu(3, 3, bm3) == 0);
    REQUIRE(k78557Seq.value_mod_mu(797827, 3, bm3) == 0);
    REQUIRE(k78557Seq.value_mod_mu(877532, 3, bm3) == 1);
}

TEST_CASE("Basic tests for Riesel sequence 3")
{
    auto k271129Seq = ctd::riesel{271129}.math_sequence();

    auto bm5 = ctd::compute_barrett_mu<ctd::mu_both_t>(5);
    REQUIRE(k271129Seq.value_mod_mu(0, 5, bm5) == 271128 % 5);
    REQUIRE(k271129Seq.value_mod_mu(1, 5, bm5) == 2);
    REQUIRE(k271129Seq.value_mod_mu(2, 5, bm5) == 0);
    REQUIRE(k271129Seq.value_mod_mu(3, 5, bm5) == 1);
    REQUIRE(k271129Seq.value_mod_mu(4, 5, bm5) == 3);
    REQUIRE(k271129Seq.value_mod_mu(6700224, 5, bm5) == 3);
    REQUIRE(k271129Seq.value_mod_mu(4924567891ull, 5, bm5) == 1);
}

TEST_CASE("Modular recurrence tests for Riesel sequence 1")
{
    auto k271129Seq = ctd::riesel{271129}.math_sequence();

    auto bm5 = ctd::compute_barrett_mu<ctd::mu_both_t>(5);
    auto n0 = 6700224u;
    auto r0 = k271129Seq.value_mod_mu(n0, 5, bm5);
    REQUIRE(k271129Seq.next_value_mod(r0, n0, 5) == 2);
    REQUIRE(k271129Seq.next_value_mod(2, n0 + 1, 5) == 0);
    REQUIRE(k271129Seq.next_value_mod(0, n0 + 2, 5) == 1);
    REQUIRE(k271129Seq.next_value_mod(1, n0 + 3, 5) == 3);
    REQUIRE(k271129Seq.next_value_mod(3, n0 + 4, 5) == 2);
    REQUIRE(k271129Seq.next_value_mod(2, n0 + 5, 5) == 0);
    REQUIRE(k271129Seq.next_value_mod(0, n0 + 6, 5) == 1);
    REQUIRE(k271129Seq.next_value_mod(1, n0 + 7, 5) == 3);
}
