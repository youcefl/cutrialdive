/*
* Created on 2026.06.03
* Copyright (c) Youcef Lemsafer
*/
#include <sstream>
#include <catch2/catch_test_macros.hpp>

#include "trial_factoring.hpp"
#include "tf.hpp"


using namespace cutrialdive;

namespace cutrialdive::tests {

    template <typename PeanoTFunc>
    void test_peano(PeanoTFunc tfPeano)
    {
        constexpr uint64_t n0 = 2956708207;
        constexpr uint64_t n1 = 2956708707;
        constexpr auto maxFactorsPerNumber = 1024;
        constexpr auto isProgressEnabled = false;
        trial_factoring_options tfOpts{n0, n1, 0, 1u << 24, {},
                maxFactorsPerNumber, isProgressEnabled};
        std::ostringstream ostr;
        auto results = tfPeano(tfOpts, ostr);

        REQUIRE(results.size() == n1 - n0);
        auto factors = results[0];
        REQUIRE(factors.size() == 0); // 2956708207 is prime
        factors = results[1];
        REQUIRE(factors.size() == 3); // 2956708208 = 2^4*43*4297541
        REQUIRE(factors[0].prime == 2);
        REQUIRE(factors[1].prime == 43);
        REQUIRE(factors[2].prime == 4297541);
        factors = results[2];
        REQUIRE(factors.size() == 6); // 2956708209 = 3*7*13*107*127*797
        REQUIRE(factors[0].prime == 3);
        REQUIRE(factors[1].prime == 7);
        REQUIRE(factors[2].prime == 13);
        REQUIRE(factors[3].prime == 107);
        REQUIRE(factors[4].prime == 127);
        REQUIRE(factors[5].prime == 797);
        factors = results[n1 - 1 - n0];
        REQUIRE(factors.size() == 3); // 2956708706 = 2*7^3*4310071
        REQUIRE(factors[0].prime == 2);
        REQUIRE(factors[1].prime == 7);
        REQUIRE(factors[2].prime == 4310071);
    }

}

TEST_CASE("TF Peano basic")
{
    tests::test_peano(tests::basic_tf_peano);
}

TEST_CASE("TF Peano with val, valmod_2, valmod_mu64")
{
    tests::test_peano(tests::value_mod_2_mu64_tf_peano);
}

TEST_CASE("TF Peano with val, valmod_2, valmod_mu128, next_val_mod_mu128")
{
    tests::test_peano(tests::value_mod_2_mu128_tf_peano);
}


TEST_CASE("TF Mersenne")
{
    constexpr auto maxFactorsPerNumber = 1024;
    constexpr auto isProgressEnabled = true;
    trial_factoring_options tfOpts{1, 101, 0, 1u << 16, {},
            maxFactorsPerNumber, isProgressEnabled};
    std::ostringstream ostr;

    auto results = trial_factor(num_seq_spec{num_seq_id::mersenne}, tfOpts, ostr);

    REQUIRE(results.n0() == 1);
    REQUIRE(results.size() == 100);
    auto m97Factors = results[97 - 1];
    REQUIRE(m97Factors.size() == 1);
    REQUIRE(m97Factors[0].prime == 11447);
    REQUIRE(m97Factors[0].exponent == 1);
}


