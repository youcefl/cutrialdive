/*
* MIT License
* Created on 2026.07.13
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#include <catch2/catch_test_macros.hpp>
#include "number_sequences/lucas.hpp"
#include "barrett_reciprocals.hpp"

namespace ctd = cutrialdive;


TEST_CASE("Basic tests for Lucas numbers 1")
{
    auto luc = ctd::lucas{};

    REQUIRE(luc.value(0) == 2);
    REQUIRE(luc.value(1) == 1);
    REQUIRE(luc.value(2) == 3);
    REQUIRE(luc.value(3) == 4);
    REQUIRE(luc.value(4) == 7);
    REQUIRE(luc.value(5) == 11);
    REQUIRE(luc.value(6) == 18);
    REQUIRE(luc.value(7) == 29);
    REQUIRE(luc.value(8) == 47);
    REQUIRE(luc.value(9) == 76);
    REQUIRE(luc.value(10) == 123);
    REQUIRE(luc.value(11) == 199);
}

TEST_CASE("Basic tests for Lucas numbers 2")
{
    auto luc = ctd::lucas{};

    REQUIRE(luc.value(36) == 33385282);
    REQUIRE(luc.value(37) == 54018521);
    REQUIRE(luc.value(38) == 87403803);
    REQUIRE(luc.value(69) == 263115950957276);
}

TEST_CASE("Lucas modular tests 1")
{
    auto luc = ctd::lucas{};

    REQUIRE(luc.value_mod_2(0) == 0);
    REQUIRE(luc.value_mod_2(1) == 1);
    REQUIRE(luc.value_mod_2(2) == 1);
    REQUIRE(luc.value_mod_2(3) == 0);
    REQUIRE(luc.value_mod_2(4) == 1);
    REQUIRE(luc.value_mod_2(5) == 1);
    REQUIRE(luc.value_mod_2(6) == 0);
    REQUIRE(luc.value_mod_2(7) == 1);
}

TEST_CASE("Lucas modular tests 2")
{
    auto luc = ctd::lucas{};

    REQUIRE(luc.next_value_mod_2(1, 28) == 1);
    REQUIRE(luc.next_value_mod_2(1, 29) == 0);
    REQUIRE(luc.next_value_mod_2(0, 30) == 1);
    REQUIRE(luc.next_value_mod_2(1, 70) == 1);
    REQUIRE(luc.next_value_mod_2(1, 71) == 0);
}

TEST_CASE("Lucas modular tests 3")
{
    auto luc = ctd::lucas{};
    {
        auto st = luc.make_state(1, 11, ctd::compute_barrett_mu_numseq<ctd::lucas>(11));

        REQUIRE(st.previous == 2);
        REQUIRE(st.current == 1);
    }
    {
        auto st = luc.make_state(11, 11, ctd::compute_barrett_mu_numseq<ctd::lucas>(11));

        REQUIRE(st.previous == 2);
        REQUIRE(st.current == 1);
    }
    {
        auto st = luc.make_state(64, 13, ctd::compute_barrett_mu_numseq<ctd::lucas>(13));

        REQUIRE(st.previous == 3);
        REQUIRE(st.current == 8);
    }
}

TEST_CASE("Lucas modular tests 4")
{
    auto luc = ctd::lucas{};
    auto st = luc.make_state(104, 97, ctd::compute_barrett_mu_numseq<ctd::lucas>(97));

    REQUIRE(st.previous == 86);
    REQUIRE(st.current == 79);
    REQUIRE(luc.value_mod(104, 97, st) == 79);
    REQUIRE(luc.next_value_mod(79, 104, 97, st) == 68);
}

