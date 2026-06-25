/*
* MIT License
* Created on 2026.05.23
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#include <array>

#include <catch2/catch_test_macros.hpp>
#include "modular_arithmetic_detail.hpp"

using namespace cutrialdive;

namespace {
    struct mu64_case {
        uint64_t p;
        uint64_t expected_mu;
    };
    struct mu128_case {
        uint64_t p;
        __uint128_t expected_mu;
    };
    struct mul128hi_case {
        __uint128_t a, b, expected_mul128hi;
    };
    template <typename T>
    inline
    T high64(__uint128_t x) { return T((x) >> 64); }
    template <typename T>
    inline
    T low64(__uint128_t x) { return T((x) & ~T{}); };
}


TEST_CASE("basic arithmetic")
{
    // I hope this one passes
    REQUIRE(2 + 2 == 4);
}

TEST_CASE("modular propagation sanity")
{
    constexpr uint64_t p = 101;

    uint64_t r = 1;

    for (uint64_t i = 0; i < 1000; ++i)
    {
        r = (r * 10) % p;
    }

    REQUIRE(r < p);
}

TEST_CASE("modular_arithmetic", "[barrett][mu64]")
{
    // Test cases generated using gen_mu64_testcases()
    auto const cases = std::to_array<mu64_case>({
        { 3ull, 0x5555555555555555ull },
        { 5ull, 0x3333333333333333ull },
        { 7ull, 0x2492492492492492ull },
        { 9ull, 0x1c71c71c71c71c71ull },
        { 15ull, 0x1111111111111111ull },
        { 17ull, 0xf0f0f0f0f0f0f0full },
        { 19ull, 0xd79435e50d79435ull },
        { 121ull, 0x21d9ead7cd391fbull },
        { 257ull, 0xff00ff00ff00ffull },
        { 521ull, 0x7dc9f3397d4c29ull },
        { 991ull, 0x4221950db0f3dbull },
        { 65535ull, 0x1000100010001ull },
        { 65537ull, 0xffff0000ffffull },
        { 988027ull, 0x10fb039bc7c1ull },
        { 4294967291ull, 0x100000005ull },
        { 18446744030759878681ull, 0x1ull },
        { 18446744065119617025ull, 0x1ull },
        { 18446744073709551557ull, 0x1ull }
    });

    for(auto const & tc : cases)
    {
        CAPTURE(tc.p);

        REQUIRE(mu64(tc.p) == tc.expected_mu);
    }        
}

TEST_CASE("modular_arithmetic", "[barrett][mu128]")
{
    // Test cases generated using gen_mu128_testcases()
    auto const cases = std::to_array<mu128_case>({
        { 3ull, (__uint128_t{0x5555555555555555} << 64) | 0x5555555555555555ull },
        { 5ull, (__uint128_t{0x3333333333333333} << 64) | 0x3333333333333333ull },
        { 7ull, (__uint128_t{0x2492492492492492} << 64) | 0x4924924924924924ull },
        { 9ull, (__uint128_t{0x1c71c71c71c71c71} << 64) | 0xc71c71c71c71c71cull },
        { 15ull, (__uint128_t{0x1111111111111111} << 64) | 0x1111111111111111ull },
        { 17ull, (__uint128_t{0xf0f0f0f0f0f0f0f} << 64) | 0xf0f0f0f0f0f0f0full },
        { 19ull, (__uint128_t{0xd79435e50d79435} << 64) | 0xe50d79435e50d794ull },
        { 121ull, (__uint128_t{0x21d9ead7cd391fb} << 64) | 0xc4c2a50658dc0876ull },
        { 255ull, (__uint128_t{0x101010101010101} << 64) | 0x101010101010101ull },
        { 257ull, (__uint128_t{0xff00ff00ff00ff} << 64) | 0xff00ff00ff00ffull },
        { 521ull, (__uint128_t{0x7dc9f3397d4c29} << 64) | 0x4643cedd1cfd8b0eull },
        { 991ull, (__uint128_t{0x4221950db0f3db} << 64) | 0xd5a27c833aa3c72bull },
        { 65535ull, (__uint128_t{0x1000100010001} << 64) | 0x1000100010001ull },
        { 65537ull, (__uint128_t{0xffff0000ffff} << 64) | 0xffff0000ffffull },
        { 988027ull, (__uint128_t{0x10fb039bc7c1} << 64) | 0x4fd0259d04ae4129ull },
        { 4294967291ull, (__uint128_t{0x100000005} << 64) | 0x190000007dull },
        { 4294967295ull, (__uint128_t{0x100000001} << 64) | 0x100000001ull },
        { 4294967297ull, (__uint128_t{0xffffffff} << 64) | 0xffffffffull },
        { 18446744030759878681ull, (__uint128_t{0x1} << 64) | 0xa0000004bull },
        { 18446744065119617025ull, (__uint128_t{0x1} << 64) | 0x200000003ull },
        { 18446744073709551557ull, (__uint128_t{0x1} << 64) | 0x3bull },
        { 18446744073709551615ull, (__uint128_t{0x1} << 64) | 0x1ull }
    });

    for(auto const & tc : cases)
    {
        CAPTURE(tc.p);

        REQUIRE(mu128(tc.p) == tc.expected_mu);
    }        
}

TEST_CASE("modular_arithmetic", "[mulhi128]")
{
    // Test cases generated using gen_mulhi128_testcases()
    auto const cases = std::to_array<mul128hi_case>({
{ 0x0ull, 0x1ull, 0x0ull },
{ 0x1ull, 0x0ull, 0x0ull },
{ 0x2ull, 0x4ull, 0x0ull },
{ 0x2ull, 0x9ull, 0x0ull },
{ 0x4ull, 0x2ull, 0x0ull },
{ 0x9ull, 0x2ull, 0x0ull },
{ 0xfffull, 0x1000000ull, 0x0ull },
{ 0x7ffffffffffull, 0x80000000001ull, 0x0ull },
{ 0xffffffffffffffffull, 0xffffffffffffffffull, 0x0ull },
{ 0xffffffffffffffffull, (__uint128_t{0x1ull} << 64) | 0x0ull, 0x0ull },
{ (__uint128_t{0x1ull} << 64) | 0x0ull, 0xffffffffffffffffull, 0x0ull },
{ (__uint128_t{0x1ull} << 64) | 0x0ull, (__uint128_t{0x1ull} << 64) | 0x0ull, 0x1ull },
{ (__uint128_t{0x1ull} << 64) | 0xfffffffe00000001ull, (__uint128_t{0x1ull} << 64) | 0xffffffffffffffffull, 0x3ull },
{ 0x20000000000000ull, (__uint128_t{0x800ull} << 64) | 0x0ull, 0x1ull },
{ (__uint128_t{0x2000ull} << 64) | 0x0ull, 0x8000000000000ull, 0x1ull },
{ (__uint128_t{0x840000000ull} << 64) | 0x891087b8b0347115ull, (__uint128_t{0x2000000000ull} << 64) | 0x1fe0000001ull,
    (__uint128_t{0x108ull} << 64) | 0x112210f81cull },
{ (__uint128_t{0xeb7fb06849c5843ull} << 64) | 0x71dcbfb976e9a33ull, (__uint128_t{0x5cf3a89821482e95ull} << 64) | 0xe382ea861ae682b2ull,
    (__uint128_t{0x558208b1abb40cdull} << 64) | 0x53213f87b19e3cb6ull },
{ (__uint128_t{0x1fffffffffffffffull} << 64) | 0xffffffffffffffffull, (__uint128_t{0x8ull} << 64) | 0x0ull, 0xffffffffffffffffull },
{ (__uint128_t{0x984f9178540196a9ull} << 64) | 0xf64fdb69582b2964ull, (__uint128_t{0xbc529db2d7a28b58ull} << 64) | 0xda7b77de9bb5098dull,
    (__uint128_t{0x700b96242a118cc4ull} << 64) | 0xa3da95e27ef5ccfcull },
{ (__uint128_t{0xabb62948a7fb7b3full} << 64) | 0x8e6cb2eb530a1cdaull, (__uint128_t{0x9b6ad2d7834f3926ull} << 64) | 0xb156d7c99f3a58d7ull,
    (__uint128_t{0x683ef1d811b8723eull} << 64) | 0x516b75e32917d8faull },
{ (__uint128_t{0x668eea7b862b221full} << 64) | 0x82813d5aa020b600ull, (__uint128_t{0xf99246d69d14700eull} << 64) | 0x116bf17a60eed1c9ull,
    (__uint128_t{0x63fb9bf4ee761898ull} << 64) | 0x4b9bcf348c418f13ull },
{ (__uint128_t{0xf0370ab5d6a70b04ull} << 64) | 0x92d71d59775f32bbull, (__uint128_t{0x597c74ada5452af8ull} << 64) | 0x3c4b13226d6de952ull,
    (__uint128_t{0x53f7eade48d50ef2ull} << 64) | 0xb928811d7061773full },
{ (__uint128_t{0x613f90a3232a3524ull} << 64) | 0x4391fb0150ddc03full, (__uint128_t{0xc7852de57fde44ecull} << 64) | 0x98cd1766b758b5aaull,
    (__uint128_t{0x4bcb00e44ff31b68ull} << 64) | 0xcb713fa4a6bb6f4cull },
{ (__uint128_t{0xfdafdafdafdafdafull} << 64) | 0xdafdafdafdafdafdull, (__uint128_t{0xafdafdafdafdafdaull} << 64) | 0xfdafdafdafdafdafull,
    (__uint128_t{0xae4439d8f6e503a9ull} << 64) | 0x902e5c3b590ee682ull },
{ (__uint128_t{0xfdfdfdfdfdfdfdfdull} << 64) | 0xfdfdfdfdfdfdfdfdull, (__uint128_t{0xfdfdfdfdfdfdfdfdull} << 64) | 0xfdfdfdfdfdfdfdfdull,
    (__uint128_t{0xfc0004080c101418ull} << 64) | 0x1c2024282c303436ull },
{ (__uint128_t{0xffffffffffffffffull} << 64) | 0xffffffffffffff61ull, (__uint128_t{0xffffffffffffffffull} << 64) | 0xffffffffffffffffull,
    (__uint128_t{0xffffffffffffffffull} << 64) | 0xffffffffffffff60ull },
{ (__uint128_t{0xffffffffffffffffull} << 64) | 0xffffffffffffffffull, (__uint128_t{0xffffffffffffffffull} << 64) | 0xffffffffffffff61ull,
    (__uint128_t{0xffffffffffffffffull} << 64) | 0xffffffffffffff60ull },
{ (__uint128_t{0xffffffffffffffffull} << 64) | 0xffffffffffffff61ull, (__uint128_t{0xffffffffffffffffull} << 64) | 0xffffffffffffff61ull,
    (__uint128_t{0xffffffffffffffffull} << 64) | 0xfffffffffffffec2ull },
{ (__uint128_t{0xffffffffffffffffull} << 64) | 0xffffffffffffffffull, (__uint128_t{0xffffffffffffffffull} << 64) | 0xffffffffffffffffull,
    (__uint128_t{0xffffffffffffffffull} << 64) | 0xfffffffffffffffeull },
    });

    for(auto const & tc : cases)
    {
        auto hi = high64<uint64_t>;
        auto lo = low64<uint64_t>;

        CAPTURE(hi(tc.a), lo(tc.a));
        CAPTURE(hi(tc.b), lo(tc.b));
        auto actual = mul128hi(tc.a, tc.b);
        auto expected = tc.expected_mul128hi;
        CAPTURE(hi(actual), lo(actual));
        CAPTURE(hi(expected), lo(expected));

        REQUIRE(actual == expected);
    }        

}
