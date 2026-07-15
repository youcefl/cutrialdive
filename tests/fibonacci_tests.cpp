/*
* MIT License
* Created on 2026.07.13
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#include <catch2/catch_test_macros.hpp>
#include "number_sequences/fibonacci.hpp"
#include "barrett_reciprocals.hpp"

namespace ctd = cutrialdive;

TEST_CASE("Basic tests for Fibonacci numbers 1")
{
    auto fib = ctd::fibonacci{};

    REQUIRE(fib.value(0) == HgInt{});
    REQUIRE(fib.value(1) == 1);
    REQUIRE(fib.value(2) == 1);
    REQUIRE(fib.value(3) == 2);
    REQUIRE(fib.value(4) == 3);
    REQUIRE(fib.value(5) == 5);
    REQUIRE(fib.value(6) == 8);
    REQUIRE(fib.value(7) == 13);
    REQUIRE(fib.value(8) == 21);
    REQUIRE(fib.value(9) == 34);
    REQUIRE(fib.value(10) == 55);
}

TEST_CASE("Basic tests for Fibonacci numbers 2")
{
    auto fib = ctd::fibonacci{};

    REQUIRE(fib.value(29) == 514229);
    REQUIRE(fib.value(30) == 832040);
    REQUIRE(fib.value(31) == 1346269);
    REQUIRE(fib.value(82) == 61305790721611591);
    REQUIRE(fib.value(83) == 99194853094755497);
}

TEST_CASE("Fibonacci modular tests 1")
{
    auto fib = ctd::fibonacci{};

    REQUIRE(fib.value_mod_2(0) == 0);
    REQUIRE(fib.value_mod_2(1) == 1);
    REQUIRE(fib.value_mod_2(2) == 1);
    REQUIRE(fib.value_mod_2(3) == 0);
    REQUIRE(fib.value_mod_2(4) == 1);
    REQUIRE(fib.value_mod_2(5) == 1);
    REQUIRE(fib.value_mod_2(6) == 0);
    REQUIRE(fib.value_mod_2(7) == 1);
}

TEST_CASE("Fibonacci modular tests 2")
{
    auto fib = ctd::fibonacci{};

    REQUIRE(fib.value_mod_2(29) == 1);
    REQUIRE(fib.value_mod_2(30) == 0);
    REQUIRE(fib.value_mod_2(31) == 1);
    REQUIRE(fib.value_mod_2(32) == 1);
    REQUIRE(fib.value_mod_2(33) == 0);
}

TEST_CASE("Fibonacci modular tests 3")
{
    auto fib = ctd::fibonacci{};

    REQUIRE(fib.next_value_mod_2(1, 28) == 1);
    REQUIRE(fib.next_value_mod_2(1, 29) == 0);
    REQUIRE(fib.next_value_mod_2(0, 30) == 1);
    REQUIRE(fib.next_value_mod_2(1, 70) == 1);
    REQUIRE(fib.next_value_mod_2(1, 71) == 0);
}

TEST_CASE("Fibonacci modular tests 4")
{
    auto fib = ctd::fibonacci{};

    REQUIRE(fib.next_value_mod_2(1, 28) == 1);
    REQUIRE(fib.next_value_mod_2(1, 29) == 0);
    REQUIRE(fib.next_value_mod_2(0, 30) == 1);
    REQUIRE(fib.next_value_mod_2(1, 70) == 1);
    REQUIRE(fib.next_value_mod_2(1, 71) == 0);
}

TEST_CASE("Fibonacci modular tests 5")
{
    auto fib = ctd::fibonacci{};
    {
        auto st = fib.make_state(1, 11, ctd::compute_barrett_mu_numseq<ctd::fibonacci>(11));

        REQUIRE(st.previous == 0);
        REQUIRE(st.current == 1);
    }
    {
        auto st = fib.make_state(11, 11, ctd::compute_barrett_mu_numseq<ctd::fibonacci>(11));

        REQUIRE(st.previous == 0);
        REQUIRE(st.current == 1);
    }
    {
        auto st = fib.make_state(64, 13, ctd::compute_barrett_mu_numseq<ctd::fibonacci>(13));

        REQUIRE(st.previous == 0);
        REQUIRE(st.current == 8);
    }
}

TEST_CASE("Fibonacci modular tests 6")
{
    auto fib = ctd::fibonacci{};

    auto st = fib.make_state(40, 89, ctd::compute_barrett_mu_numseq<ctd::fibonacci>(89));

    REQUIRE(st.previous == 5);
    REQUIRE(st.current == 86);
    REQUIRE(fib.value_mod(40, 89, st) == 86);
    REQUIRE(fib.next_value_mod(86, 40, 89, st) == 2);
    REQUIRE(fib.next_value_mod(2, 41, 89, st) == 88);
    REQUIRE(fib.next_value_mod(88, 42, 89, st) == 1);
    REQUIRE(fib.next_value_mod(1, 43, 89, st) == 0);
    REQUIRE(fib.next_value_mod(0, 44, 89, st) == 1);
}
