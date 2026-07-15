/*
* MIT License
* Created on 2026.07.14
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#pragma once

#include <cstdint>
#include <iostream>
#include "common_defines.h"
#include "barrett_mu_types.hpp"
#include "modular_arithmetic_detail.hpp"
#include "hgint.hpp"
#include "mod_arith.hpp"

namespace cutrialdive {

    class tribonacci
    {
    public:
        struct state
        {
            uint64_t t_n;  // T(n)
            uint64_t t_n1;  // T(n+1)
            uint64_t t_n2;  // T(n+2)
        };

        using index_type = uint64_t;
        using residue_type = uint64_t;
        using mu_type = mu_both_t;
        static const bool has_state = true;
        using state_type = state;

        CUTRIALDIVE_DEVICE_AND_HOST state_type make_state(index_type n, residue_type d, mu_both_t mu);
        CUTRIALDIVE_DEVICE_AND_HOST residue_type value_mod_2(index_type n);
        CUTRIALDIVE_DEVICE_AND_HOST residue_type next_value_mod_2(residue_type r, index_type n);
        CUTRIALDIVE_DEVICE_AND_HOST residue_type value_mod(index_type n, residue_type d, state_type & st);
        CUTRIALDIVE_DEVICE_AND_HOST residue_type next_value_mod(residue_type r, index_type n, residue_type d, state_type & st);

        HgInt value(index_type n);
        std::ostream& print_value(index_type n, std::ostream & out);
        static char const * short_name();
    };

}

namespace cutrialdive {

    // From Alexander Samokrutov (June 2026, see formula section of https://oeis.org/A000073):
    // a(n+m) = (a(m+2)-a(m+1)-a(m))*a(n) + a(m)*a(n+2) + (a(m+1)-a(m))*a(n+1)
    // a(n+m+1) = a(m)*a(n) + a(m+1)*a(n+2) - (a(m+1)-a(m+2))*a(n+1)
    // a(n+m+2) = a(m+1)*a(n) + a(m+2)*a(n+2) + (a(m+1)+a(m))*a(n+1).

    // setting n = m we get

    // T_2n = (T_{n+2} - T_{n+1} - T_n)*T_n + T_n*T_{n+2} + (T_{n+1} - T_n)*T_{n+1}
    //     = 2T_{n+2}*T_n + T_{n+1}^2 - 2T_{n+1}*T_n - T_n^2
    // T_{2n+1} = T_n^2 + T_{n+1}*T_{n+2} + (T_{n+2} - T_{n+1})*T_{n+1}
    //         = T_n^2 + 2*T_{n+2}*T_{n+1} - T_{n+1}^2
    // T_{2n+2} = 2*T_{n+1}*T_n + T_{n+2}^2 + T_{n+1}^2

    // We can then use binary Lucas chain to evaluate T_n (mod d).
    template <typename ArithT, typename IndexT>
    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    void tribonacci_at(ArithT arith,
        IndexT n,
        typename ArithT::value_type & vn,
        typename ArithT::value_type & vn1,
        typename ArithT::value_type & vn2
    )
    {
        using V = typename ArithT::value_type;
        vn = V{};
        vn1 = V{};
        vn2 = V{1};
        if(!n) {
            return;
        }
        for(IndexT m = IndexT{1} << (bit_length(n) - 1); m; m >>= 1) {
            auto sqr_vn = arith.sqr(vn);
            auto sqr_vn1 = arith.sqr(vn1);
            auto two_vn1_vn = arith.mul(arith.mul_by_2(vn1), vn);
            auto tvn = arith.sub(arith.sub(arith.add(
                                    arith.mul(arith.mul_by_2(vn2), vn),
                                    sqr_vn1),
                                two_vn1_vn),
                                sqr_vn);
            auto tvn1 = arith.sub(arith.add(sqr_vn, 
                                    arith.mul(arith.mul_by_2(vn2), vn1)),
                                sqr_vn1);
            vn2 = arith.add(arith.add(two_vn1_vn, arith.sqr(vn2)), sqr_vn1);
            vn = std::move(tvn);
            vn1 = std::move(tvn1);
            if(m & n) {
                auto vn3 = arith.add(arith.add(vn, vn1), vn2);
                vn = std::move(vn1);
                vn1 = std::move(vn2);
                vn2 = std::move(vn3);
            }
        }
    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    tribonacci::residue_type tribonacci::value_mod_2(index_type n)
    {
        return (n & 3) >> 1;
    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    tribonacci::residue_type tribonacci::next_value_mod_2(residue_type r, index_type n)
    {
        (void)r;
        return value_mod_2(n + 1);
    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    tribonacci::state_type tribonacci::make_state(index_type n, residue_type d, mu_both_t mu)
    {
        state_type st;
        tribonacci_at(mod_arith{d, mu}, n, st.t_n, st.t_n1, st.t_n2);
        return st;
    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    tribonacci::residue_type tribonacci::value_mod(index_type n, residue_type d, state_type & st)
    {
        (void)n;
        (void)d;
        return st.t_n;
    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    tribonacci::residue_type tribonacci::next_value_mod(residue_type r, index_type n, residue_type d, state_type & st)
    {
        uint64_t t_n3 = addmod(addmod(st.t_n, st.t_n1, d), st.t_n2, d);
        st.t_n = st.t_n1;
        st.t_n1 = st.t_n2;
        st.t_n2 = t_n3;
        return st.t_n;
    }

}
