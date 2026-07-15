/*
* MIT License
* Created on 2026.07.15
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#include <catch2/catch_test_macros.hpp>
#include "number_sequences/tribonacci.hpp"
#include "barrett_reciprocals.hpp"

namespace ctd = cutrialdive;

// See https://oeis.org/A000073 for sequence definition and values

TEST_CASE("Basic tests for tribonacci numbers 1")
{
    auto trib = ctd::tribonacci{};

    REQUIRE(trib.value(0) == HgInt{});
    REQUIRE(trib.value(1) == HgInt{});
    REQUIRE(trib.value(2) == 1);
    REQUIRE(trib.value(3) == 1);
    REQUIRE(trib.value(4) == 2);
    REQUIRE(trib.value(5) == 4);
    REQUIRE(trib.value(6) == 7);
    REQUIRE(trib.value(7) == 13);
    REQUIRE(trib.value(8) == 24);
    REQUIRE(trib.value(9) == 44);
    REQUIRE(trib.value(10) == 81);
}

TEST_CASE("Basic tests for tribonacci numbers 2")
{
    auto trib = ctd::tribonacci{};

    REQUIRE(trib.value(52) == 10562230626642);
    REQUIRE(trib.value(53) == 19426970897100);
    REQUIRE(trib.value(54) == 35731770264967);
    REQUIRE(trib.value(55) == 65720971788709);
}

TEST_CASE("Basic tests for tribonacci numbers 3")
{
    auto trib = ctd::tribonacci{};

    REQUIRE(trib.value(199) == HgInt{"8457146198099353452163566448904305495553772125703149"});
    REQUIRE(trib.value(200) == HgInt{"15555116989073938986569525465884451018665640926743832"});
    REQUIRE(trib.value(201) == HgInt{"28610320653810477165032088685001500201865067503083660"});
}

TEST_CASE("Tribonacci modular tests 1")
{
    auto trib = ctd::tribonacci{};

    REQUIRE(trib.value_mod_2(0) == 0);
    REQUIRE(trib.value_mod_2(1) == 0);
    REQUIRE(trib.value_mod_2(2) == 1);
    REQUIRE(trib.value_mod_2(3) == 1);
    REQUIRE(trib.value_mod_2(4) == 0);
    REQUIRE(trib.value_mod_2(5) == 0);
    REQUIRE(trib.value_mod_2(6) == 1);
    REQUIRE(trib.value_mod_2(7) == 1);
}

TEST_CASE("Tribonacci modular tests 2")
{
    auto trib = ctd::tribonacci{};

    REQUIRE(trib.value_mod_2(199) == 1);
    REQUIRE(trib.value_mod_2(200) == 0);
    REQUIRE(trib.value_mod_2(201) == 0);
    REQUIRE(trib.value_mod_2(1000) == 0);
    REQUIRE(trib.value_mod_2(1001) == 0);
    REQUIRE(trib.value_mod_2(1002) == 1);
    REQUIRE(trib.value_mod_2(1003) == 1);
    REQUIRE(trib.value_mod_2(1004) == 0);
    REQUIRE(trib.value_mod_2(1005) == 0);
    REQUIRE(trib.value_mod_2(1006) == 1);
    REQUIRE(trib.value_mod_2(1007) == 1);
}

TEST_CASE("Tribonacci modular tests 3")
{
    auto trib = ctd::tribonacci{};

    REQUIRE(trib.next_value_mod_2(1, 199) == 0);
    REQUIRE(trib.next_value_mod_2(0, 200) == 0);
    REQUIRE(trib.next_value_mod_2(0, 201) == 1);
}

TEST_CASE("Tribonacci modular tests 5")
{
    auto trib = ctd::tribonacci{};
    {
        auto st = trib.make_state(1, 11, ctd::compute_barrett_mu_numseq<ctd::tribonacci>(11));

        REQUIRE(st.t_n == 0);
        REQUIRE(st.t_n1 == 1);
        REQUIRE(st.t_n2 == 1);
    }
    {
        auto st = trib.make_state(7, 11, ctd::compute_barrett_mu_numseq<ctd::tribonacci>(11));

        REQUIRE(st.t_n == 2);
        REQUIRE(st.t_n1 == 2);
        REQUIRE(st.t_n2 == 0);
    }
    {
        uint64_t d = (1u << 31) - 1;
        auto st = trib.make_state(64, d, ctd::compute_barrett_mu_numseq<ctd::tribonacci>(d));

        REQUIRE(st.t_n == 15832480722303616 % d);
        REQUIRE(st.t_n1 == 29120472094716576 % d);
        REQUIRE(st.t_n2 == 53560898629395777 % d);
    }
}
