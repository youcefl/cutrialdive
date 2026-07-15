/*
* MIT License
* Created on 2026.07.11
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

namespace cutrialdive {

    class fibonacci
    {
    public:
        struct state
        {
            uint64_t previous;
            uint64_t current;
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

        HgInt value(index_type n, HgInt & nextValue);
    };

}

namespace cutrialdive {

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    fibonacci::residue_type fibonacci::value_mod_2(index_type n)
    {
        return n % 3 == 0 ? 0 : 1;

    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    fibonacci::residue_type fibonacci::next_value_mod_2(residue_type r, index_type n)
    {
        (void)r;
        return (n + 1) % 3 == 0 ? 0 : 1;
    }


    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    fibonacci::state_type   fibonacci::make_state(index_type n, residue_type d, mu_both_t mu)
    {
        state_type st{};
        if(n <= 1) {
            st.current = n;
            return st;
        }
        n -= 1;
        residue_type f0 = 0, f1 = 1;
        for(index_type m = index_type{1} << (bit_length(n) - 1); m; m >>= 1) {
            if(m & n) {
                auto t0 = mulmod(submod(doublemod(f1, d), f0, d), f0, d, mu.mu64, mu.mu128);
                f0 = addmod(mulmod(f1, f1, d, mu.mu64, mu.mu128), mulmod(f0, f0, d, mu.mu64, mu.mu128), d);
                f1 = addmod(f0, t0, d);
            } else {
                auto t0 = f0;
                f0 = mulmod(submod(doublemod(f1, d), f0, d), f0, d, mu.mu64, mu.mu128);
                f1 = addmod(mulmod(f1, f1, d, mu.mu64, mu.mu128), mulmod(t0, t0, d, mu.mu64, mu.mu128), d);
            }
        }
        st.previous = f0;
        st.current = f1;
        return st;
    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    fibonacci::residue_type fibonacci::value_mod(index_type n, residue_type d, state_type & st)
    {
        (void)n;
        (void)d;
        return st.current;
    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    fibonacci::residue_type fibonacci::next_value_mod(residue_type r, index_type n, residue_type d, state_type & st)
    {
        auto result = addmod(r, st.previous, d);
        st.previous = r;
        return result;
    }


}

