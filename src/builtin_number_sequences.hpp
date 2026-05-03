/*
* Creation date: 2026.04.29
* Created by Youcef Lemsafer
*/
#pragma once

#include "hgint.hpp"
#include "number_sequence.hpp"

namespace cutrialdive {

    template <>
    struct number_sequence<mode_flag::smarandache>
    {
        using index_type = uint64_t;
        using value_type = HgInt;
        using residue_type = uint64_t;
        static value_type value(index_type n);
        static residue_type next_value_mod_mu(residue_type v_n_mod_p, index_type n, residue_type mu_p);
    };

    template <>
    struct number_sequence<mode_flag::mersenne>
    {
        using index_type = uint64_t;
        using value_type = HgInt;
        using residue_type = uint64_t;
        static value_type value(index_type n);
        static residue_type next_value_mod_mu(residue_type v_n_mod_p, index_type n, residue_type mu_p);
    };

    inline uint64_t number_sequence<mode_flag::smarandache>::next_value_mod_mu(
        uint64_t v_n_mod_p,
        uint64_t n,
        uint64_t mu_p
    )
    {
        auto m = n + 1;
        uint64_t shift;
        if(m >= 10'000 && m < 100'000) {
            shift = 100'000;
        } else if(m >= 100'000 && m < 1'000'000) {
            shift = 1'000'000;
        } else if(m >= 1'000'000 && m < 10'000'000) {
            shift = 10'000'000;
        }

        return (((__uint128_t(v_n_mod_p) * shift) + m) * mu_p) >> 64;
    }

}
