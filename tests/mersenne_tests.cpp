/*
* Created on 2026.06.02
* Copyright (c) Youcef Lemsafer
*/
#include <catch2/catch_test_macros.hpp>

#include "number_sequences/mersenne.hpp"
#include "barrett_reciprocals.hpp"

using cutrialdive::mersenne;
using cutrialdive::compute_barrett_mu;
using cutrialdive::barrett_mu_both_t;


TEST_CASE("Basic Mersenne number tests 1", "[mersenne]")
{
    REQUIRE(mersenne::value(0) == mersenne::value_type{uint64_t{0}});
    REQUIRE(mersenne::value(1) == mersenne::value_type{uint64_t{1}});
    REQUIRE(mersenne::value(2) == mersenne::value_type{uint64_t{3}});
    REQUIRE(mersenne::value(3) == mersenne::value_type{uint64_t{7}});
}

TEST_CASE("Basic Mersenne number tests 2", "[mersenne]")
{
    REQUIRE(mersenne::value(127) == (mersenne::value_type{uint64_t{1}} << 127) - 1);
    REQUIRE(mersenne::value(521) == (mersenne::value_type{uint64_t{1}} << 521) - 1);
    REQUIRE(mersenne::value(991) == (mersenne::value_type{uint64_t{1}} << 991) - 1);
}

TEST_CASE("Basic Mersenne number modular tests 1", "[mersenne][modular]")
{
    REQUIRE(mersenne::value_mod_2(0) == uint64_t{0});
    REQUIRE(mersenne::value_mod_2(1) == uint64_t{1});
    REQUIRE(mersenne::value_mod_2(2) == uint64_t{1});
    REQUIRE(mersenne::value_mod_2(3) == uint64_t{1});
    REQUIRE(mersenne::value_mod_2(4) == uint64_t{1});
}

TEST_CASE("Basic Mersenne number modular tests 2", "[mersenne][modular]")
{
    REQUIRE(mersenne::next_value_mod_2(0, 0) == uint64_t{1});
    REQUIRE(mersenne::next_value_mod_2(1, 1) == uint64_t{1});
    REQUIRE(mersenne::next_value_mod_2(1, 2) == uint64_t{1});
    REQUIRE(mersenne::next_value_mod_2(1, 3) == uint64_t{1});
    REQUIRE(mersenne::next_value_mod_2(1, 4) == uint64_t{1});
    REQUIRE(mersenne::next_value_mod_2(1, 5) == uint64_t{1});
}

TEST_CASE("Mersenne number modular tests", "[mersenne][modular]")
{
    auto mu = compute_barrett_mu<barrett_mu_both_t>;
    
    REQUIRE(mersenne::value_mod_mu(21, 7, mu(7)) == 0);
    REQUIRE(mersenne::value_mod_mu(21, 127, mu(127)) == 0);
    REQUIRE(mersenne::value_mod_mu(35, 31, mu(31)) == 0);
    REQUIRE(mersenne::value_mod_mu(37, 3, mu(3)) == 1);
    REQUIRE(mersenne::value_mod_mu(257, 3, mu(3)) == 1);
    REQUIRE(mersenne::value_mod_mu(257, 5, mu(5)) == 1);
    REQUIRE(mersenne::value_mod_mu(257, 7, mu(7)) == 3);
    REQUIRE(mersenne::value_mod_mu(67, 193707721, mu(193707721)) == 0);
    REQUIRE(mersenne::value_mod_mu(521, 5, mu(5)) == 1);
    REQUIRE(mersenne::value_mod_mu(991, 5, mu(5)) == 2);
    REQUIRE(mersenne::value_mod_mu(1279, 59, mu(59)) == 7);
    REQUIRE(mersenne::value_mod_mu(131071, ~uint64_t{}, mu(~uint64_t{})) == (1ull << 63) - 1);
    REQUIRE(mersenne::value_mod_mu(131071, ~uint32_t{}, mu(~uint32_t{})) == (1ull << 31) - 1);
    REQUIRE(mersenne::value_mod_mu(61*524287, (1ull << 61) - 1, mu((1ull << 61) - 1)) == 0);
}

