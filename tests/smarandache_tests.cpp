/*
* MIT License
* Created on 2026.05.26
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#include <array>

#include <catch2/catch_test_macros.hpp>
#include "number_sequences/smarandache.hpp"
#include "barrett_reciprocals.hpp"

namespace ctd = cutrialdive;
using IntT = decltype(ctd::smarandache<10>::value(0));

namespace {
    struct sm_value_case {
        uint64_t n;
        IntT expected_sm;
    };

    /// Naive version of Smarandache base @param base
    IntT smarandache(uint64_t base, uint64_t n)
    {
        auto v = IntT{0ul};
        for(uint64_t i = 1; i <= n; ++i) {
            uint64_t exponent = 1;
            uint64_t m = i;
            while(m /= base) {
                ++exponent;
            }
            v = v * pow(IntT{base}, exponent) + i;
        }
        return v;
    }

    // Shortener for the Barrett multiplier used by the smarandache class
    // d must be odd
    auto mu(uint64_t d)
    {
        return ctd::compute_barrett_mu<ctd::mu_both_t>(d);
    }
}

TEST_CASE("basic Smarandache sequence test")
{
    auto cases = std::to_array<sm_value_case>({
        {0, IntT{uint64_t{0}}},
        {1, IntT{uint64_t{1}}},
        {2, IntT{uint64_t{12}}}
    });

    for(auto const & tc : cases)
    {
        CAPTURE(tc.n);

        REQUIRE(ctd::smarandache<10>::value(tc.n) == tc.expected_sm);
    }
}

TEST_CASE("Smarandache base 10")
{
    auto cases = std::to_array<sm_value_case>({
        {3, IntT{uint64_t{123}}},
        {4, IntT{uint64_t{1234}}},
        {5, IntT{uint64_t{12345}}},
        {6, IntT{uint64_t{123456}}},
        {7, IntT{uint64_t{1234567}}},
        {10, IntT{"12345678910"}},
        {11, IntT{"1234567891011"}},
        {12, IntT{"123456789101112"}},
        {99, smarandache(10, 99)},
        {100, smarandache(10, 100)},
        {101, smarandache(10, 101)},
        {991, smarandache(10, 991)},
        {999, smarandache(10, 999)},
        {1000, smarandache(10, 1000)},
        {1001, smarandache(10, 1001)}
    });

    for(auto const & tc : cases)
    {
        CAPTURE(tc.n);

        REQUIRE(ctd::smarandache<10>::value(tc.n) == tc.expected_sm);
    }
}

TEST_CASE("Smarandache base 2")
{
    auto cases = std::to_array<sm_value_case>({
        {0, IntT{uint64_t{0b0}}},
        {1, IntT{uint64_t{0b1}}},
        {2, IntT{uint64_t{0b110}}},
        {3, IntT{uint64_t{0b11011}}},
        {4, IntT{uint64_t{0b11011100}}},
        {5, IntT{uint64_t{0b11011100101}}},
        {6, IntT{uint64_t{0b11011100101110}}},
        {7, IntT{uint64_t{0b11011100101110111}}},
        {8, IntT{uint64_t{0b110111001011101111000}}},
        {15, smarandache(2, 15)},
        {16, smarandache(2, 16)},
        {17, smarandache(2, 17)},
        {18, smarandache(2, 18)},
        {19, smarandache(2, 19)},
        {30, smarandache(2, 30)},
        {31, smarandache(2, 31)},
        {499, smarandache(2, 499)},
        {511, smarandache(2, 511)},
        {512, smarandache(2, 512)},
        {1025, smarandache(2, 1025)},
        {3072, smarandache(2, 3072)},
        {4097, smarandache(2, 4097)},
    });

    for(auto const & tc : cases)
    {
        CAPTURE(tc.n);

        REQUIRE(ctd::smarandache<2>::value(tc.n) == tc.expected_sm);
    }
}

TEST_CASE("Smarandache base 3")
{
    auto cases = std::to_array<sm_value_case>({
        {0, IntT{uint64_t{0}}},
        {1, IntT{uint64_t{1}}},
        {2, IntT{uint64_t{3*1+2}}},
        {3, IntT{uint64_t{(3*1+2)*9+3}}},
        {4, IntT{uint64_t{((3*1+2)*9+3)*9+4}}},
        {5, IntT{uint64_t{(((3*1+2)*9+3)*9+4)*9+5}}},
        {8, smarandache(3, 8)},
        {9, smarandache(3, 9)},
        {10, smarandache(3, 10)},
        {26, smarandache(3, 26)},
        {27, smarandache(3, 27)},
        {28, smarandache(3, 28)},
        {242, smarandache(3, 242)},
        {243, smarandache(3, 243)},
        {244, smarandache(3, 244)},
    });

    for(auto const & tc : cases)
    {
        CAPTURE(tc.n);

        REQUIRE(ctd::smarandache<3>::value(tc.n) == tc.expected_sm);
    }
}    

TEST_CASE("Smarandache base 4")
{
    auto cases = std::to_array<sm_value_case>({
        {0, IntT{uint64_t{0b0}}},
        {1, IntT{uint64_t{0b1}}},
        {2, IntT{uint64_t{0b1'10}}},
        {3, IntT{uint64_t{0b1'10'11}}},
        {4, IntT{uint64_t{0b1'10'11'0100}}},
        {5, IntT{uint64_t{0b1'10'11'0100'0101}}},
        {6, IntT{uint64_t{0b1'10'11'0100'0101'0110}}},
        {15, smarandache(4, 15)},
        {16, smarandache(4, 16)},
        {17, smarandache(4, 17)},
        {63, smarandache(4, 63)},
        {64, smarandache(4, 64)},
        {65, smarandache(4, 65)},
        {255, smarandache(4, 255)},
        {256, smarandache(4, 256)},
        {257, smarandache(4, 257)},
        {2367, smarandache(4, 2367)},
        {4096, smarandache(4, 4096)},
        {4097, smarandache(4, 4097)},
    });

    for(auto const & tc : cases)
    {
        CAPTURE(tc.n);

        REQUIRE(ctd::smarandache<4>::value(tc.n) == tc.expected_sm);
    }
}    

TEST_CASE("Smarandache base 5")
{
    auto cases = std::to_array<sm_value_case>({
        {0, IntT{uint64_t{0}}},
        {1, IntT{uint64_t{1}}},
        {2, IntT{uint64_t{5+2}}},
        {3, smarandache(5, 3)},
        {4, smarandache(5, 4)},
        {5, smarandache(5, 5)},
        {6, smarandache(5, 6)},
        {124, smarandache(5, 124)},
        {125, smarandache(5, 125)},
        {126, smarandache(5, 126)},
        {624, smarandache(5, 624)},
        {625, smarandache(5, 625)},
        {626, smarandache(5, 626)},
        {2999, smarandache(5, 2999)},
        {3124, smarandache(5, 3124)},
        {3125, smarandache(5, 3125)},
        {3126, smarandache(5, 3126)},
        {3887, smarandache(5, 3887)},
    });

    for(auto const & tc : cases)
    {
        CAPTURE(tc.n);

        REQUIRE(ctd::smarandache<5>::value(tc.n) == tc.expected_sm);
    }
}

TEST_CASE("Smarandache base 8")
{
    auto cases = std::to_array<sm_value_case>({
        {0, IntT{uint64_t{00}}},
        {1, IntT{uint64_t{01}}},
        {2, IntT{uint64_t{012}}},
        {3, IntT{uint64_t{0123}}},
        {4, IntT{uint64_t{01234}}},
        {5, IntT{uint64_t{012345}}},
        {6, IntT{uint64_t{0123456}}},
        {7, IntT{uint64_t{01234567}}},
        {8, IntT{uint64_t{0123456710}}},
        {9, IntT{uint64_t{012345671011}}},
        {10, IntT{uint64_t{01234567101112}}},
        {511, smarandache(8, 511)},
        {512, smarandache(8, 512)},
        {513, smarandache(8, 513)},
        {2433, smarandache(8, 2433)},
        {4095, smarandache(8, 4095)},
        {4096, smarandache(8, 4096)},
        {4097, smarandache(8, 4097)},
    });

    for(auto const & tc : cases)
    {
        CAPTURE(tc.n);

        REQUIRE(ctd::smarandache<8>::value(tc.n) == tc.expected_sm);
    }
}

TEST_CASE("Smarandache base 16")
{
    auto cases = std::to_array<sm_value_case>({
        {0, IntT{uint64_t{0x0}}},
        {1, IntT{uint64_t{0x1}}},
        {2, IntT{uint64_t{0x12}}},
        {3, IntT{uint64_t{0x123}}},
        {4, IntT{uint64_t{0x1234}}},
        {5, IntT{uint64_t{0x12345}}},
        {6, IntT{uint64_t{0x123456}}},
        {7, IntT{uint64_t{0x1234567}}},
        {8, IntT{uint64_t{0x12345678}}},
        {9, IntT{uint64_t{0x123456789}}},
        {10, IntT{uint64_t{0x123456789a}}},
        {11, IntT{uint64_t{0x123456789ab}}},
        {12, IntT{uint64_t{0x123456789abc}}},
        {13, IntT{uint64_t{0x123456789abcd}}},
        {14, IntT{uint64_t{0x123456789abcde}}},
        {15, IntT{uint64_t{0x123456789abcdef}}},
        {255, smarandache(16, 255)},
        {256, smarandache(16, 256)},
        {257, smarandache(16, 257)},
        {3392, smarandache(16, 3392)},
        {4095, smarandache(16, 4095)},
        {4096, smarandache(16, 4096)},
        {4097, smarandache(16, 4097)},
    });

    for(auto const & tc : cases)
    {
        CAPTURE(tc.n);

        REQUIRE(ctd::smarandache<16>::value(tc.n) == tc.expected_sm);
    }
}


TEST_CASE("Smarandache base 19")
{
    auto cases = std::to_array<sm_value_case>({
        {0, IntT{uint64_t{0}}},
        {1, IntT{uint64_t{1}}},
        {2, IntT{uint64_t{19+2}}},
        {3, IntT{uint64_t{(19+2)*19+3}}},
        {4, IntT{uint64_t{((19+2)*19+3)*19+4}}},
        {5, IntT{uint64_t{(((19+2)*19+3)*19+4)*19+5}}},
        {18, smarandache(19, 18)},
        {19, smarandache(19, 19)},
        {20, smarandache(19, 20)},
        {360, smarandache(19, 360)},
        {361, smarandache(19, 361)},
        {362, smarandache(19, 362)},
        {6858, smarandache(19, 6858)},
        {6859, smarandache(19, 6859)},
        {6860, smarandache(19, 6860)},
    });

    for(auto const & tc : cases)
    {
        CAPTURE(tc.n);

        REQUIRE(ctd::smarandache<19>::value(tc.n) == tc.expected_sm);
    }
}

TEST_CASE("Smarandache expression test", "[Smarandache][base 10]")
{
    REQUIRE(ctd::smarandache_expr<10>(0) == "0");
    REQUIRE(ctd::smarandache_expr<10>(1) == "1");
    REQUIRE(ctd::smarandache_expr<10>(2) == "(10*10^2 - (10 - 1)*2 - 10)/(10 - 1)^2");
    REQUIRE(ctd::smarandache_expr<10>(3) == "(10*10^3 - (10 - 1)*3 - 10)/(10 - 1)^2");
    REQUIRE(ctd::smarandache_expr<10>(9) == "(10*10^9 - (10 - 1)*9 - 10)/(10 - 1)^2");
    REQUIRE(ctd::smarandache_expr<10>(10) ==
                "(((10*10^9 - (10 - 1)^2 - 10)/(10 - 1)^2*(10^2 - 1)^2"
                    " + 10^3 - 10 + 1)*10^2 - (10^2 - 1)*10 - 10^2)/(10^2 - 1)^2");
    REQUIRE(ctd::smarandache_expr<10>(29) ==
                "(((10*10^9 - (10 - 1)^2 - 10)/(10 - 1)^2*(10^2 - 1)^2"
                    " + 10^3 - 10 + 1)*10^40 - (10^2 - 1)*29 - 10^2)/(10^2 - 1)^2");
}

TEST_CASE("Smarandache modulo 2", "[Smarandache][base 10]")
{
    REQUIRE(ctd::smarandache<10>::value_mod_2(0) == 0);
    REQUIRE(ctd::smarandache<10>::value_mod_2(1) == 1);
    REQUIRE(ctd::smarandache<10>::value_mod_2(2) == 0);
    REQUIRE(ctd::smarandache<10>::value_mod_2(3) == 1);
    REQUIRE(ctd::smarandache<10>::value_mod_2(257) == 1);
    REQUIRE(ctd::smarandache<10>::value_mod_2(258) == 0);
    for(auto n : {258u, 640u, 1989u, 8191u})
    {
        CAPTURE(n);

        REQUIRE(
            ctd::smarandache<10>::next_value_mod_2(
                ctd::smarandache<10>::value_mod_2(n),
                n)
            ==
            ctd::smarandache<10>::value_mod_2(n + 1));
    }

    auto r = ctd::smarandache<10>::value_mod_2(0);

    for(uint64_t n = 0; n < 10000; ++n)
    {
        REQUIRE(r == ctd::smarandache<10>::value_mod_2(n));

        r = ctd::smarandache<10>::next_value_mod_2(r, n);
    }
}

TEST_CASE("Smarandache base 3 modulo 2", "[Smarandache][base 3]")
{
    REQUIRE(ctd::smarandache<3>::value_mod_2(0) == 0);
    REQUIRE(ctd::smarandache<3>::value_mod_2(1) == 1);
    REQUIRE(ctd::smarandache<3>::value_mod_2(2) == 1);
    REQUIRE(ctd::smarandache<3>::value_mod_2(3) == 0);
    for(auto n : {47u, 509u, 1024u, 6561u})
    {
        CAPTURE(n);

        REQUIRE(
            ctd::smarandache<3>::next_value_mod_2(
                ctd::smarandache<3>::value_mod_2(n),
                n)
            ==
            ctd::smarandache<3>::value_mod_2(n + 1));
    }

    auto r = ctd::smarandache<3>::value_mod_2(0);

    for(uint64_t n = 0; n < 10000; ++n)
    {
        REQUIRE(r == ctd::smarandache<3>::value_mod_2(n));

        r = ctd::smarandache<3>::next_value_mod_2(r, n);
    }
}

TEST_CASE("Smarandache base 3 modulo 3", "[Smarandache][base 3]")
{
    auto bm = mu(3);
    REQUIRE(ctd::smarandache<3>::next_value_mod_mu(0, 3, 3, bm) == 1);
    REQUIRE(ctd::smarandache<3>::next_value_mod_mu(1, 4, 3, bm) == 2);
    REQUIRE(ctd::smarandache<3>::next_value_mod_mu(2, 5, 3, bm) == 0);
    REQUIRE(ctd::smarandache<3>::next_value_mod_mu(0, 6, 3, bm) == 1);
    REQUIRE(ctd::smarandache<3>::next_value_mod_mu(1, 7, 3, bm) == 2);
    REQUIRE(ctd::smarandache<3>::next_value_mod_mu(2, 8, 3, bm) == 0);
    REQUIRE(ctd::smarandache<3>::next_value_mod_mu(0, 531441, 3, bm) == 1);
}

TEST_CASE("Smarandache base 10 modulo 5", "[Smarandache][base 10]")
{
    auto bm = mu(5);
    REQUIRE(ctd::smarandache<10>::next_value_mod_mu(0, 10, 5, bm) == 1);
    REQUIRE(ctd::smarandache<10>::next_value_mod_mu(1, 11, 5, bm) == 2);
    REQUIRE(ctd::smarandache<10>::next_value_mod_mu(2, 12, 5, bm) == 3);
    REQUIRE(ctd::smarandache<10>::next_value_mod_mu(3, 13, 5, bm) == 4);
    REQUIRE(ctd::smarandache<10>::next_value_mod_mu(4, 14, 5, bm) == 0);
    REQUIRE(ctd::smarandache<10>::next_value_mod_mu(4, 99, 5, bm) == 0);
    REQUIRE(ctd::smarandache<10>::next_value_mod_mu(0, 100, 5, bm) == 1);
}

TEST_CASE("Smarandache base 10 misc", "[Smarandache][base 10]")
{
    REQUIRE(ctd::smarandache<10>::next_value_mod_mu(13, 42, 17, mu(17)) == 0);
    REQUIRE(ctd::smarandache<10>::next_value_mod_mu(3, 42, 7, mu(7)) == 0);
    REQUIRE(ctd::smarandache<10>::next_value_mod_mu(115, 42, 7*17, mu(7*17)) == 0);
    REQUIRE(ctd::smarandache<10>::next_value_mod_mu(31, 42, 449, mu(449)) == 0);
    REQUIRE(ctd::smarandache<10>::next_value_mod_mu(2, 99, 7, mu(7)) == 0);
    REQUIRE(ctd::smarandache<10>::next_value_mod_mu(817, 99, 8171, mu(8171)) == 0);
    REQUIRE(ctd::smarandache<10>::next_value_mod_mu(0, 100, 8171, mu(8171)) == 101);
    REQUIRE(ctd::smarandache<10>::next_value_mod_mu(101, 101, 8171, mu(8171)) == (1000*101 + 102) % 8171);
    auto const m31 = (1u << 31) - 1;
    REQUIRE(ctd::smarandache<10>::next_value_mod_mu(1855760753u, 999, m31, mu(m31)) == 1201337273u);
    auto const d = ((1ull << 32) - 17)*((1ull << 32) - 5);
    REQUIRE(ctd::smarandache<10>::next_value_mod_mu(94664850192674748ull, 999, d, mu(d)) == 5864558986513650361ull);
}
