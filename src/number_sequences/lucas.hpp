/*
* MIT License
* Created on 2026.07.13
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#pragma once

#include "common_defines.h"
#include "fibonacci.hpp"


namespace cutrialdive {

    class lucas
    {
    public:
        using index_type = typename fibonacci::index_type;
        using residue_type = typename fibonacci::residue_type;
        using mu_type = typename fibonacci::mu_type;
        static const bool has_state = fibonacci::has_state;
        using state_type = typename fibonacci::state;

        CUTRIALDIVE_DEVICE_AND_HOST state_type make_state(index_type n, residue_type d, mu_both_t mu);
        CUTRIALDIVE_DEVICE_AND_HOST residue_type value_mod_2(index_type n);
        CUTRIALDIVE_DEVICE_AND_HOST residue_type next_value_mod_2(residue_type r, index_type n);
        CUTRIALDIVE_DEVICE_AND_HOST residue_type value_mod(index_type n, residue_type d, state_type & st);
        CUTRIALDIVE_DEVICE_AND_HOST residue_type next_value_mod(residue_type r, index_type n, residue_type d, state_type & st);

        HgInt value(index_type k);
        std::ostream& print_value(index_type n, std::ostream & out);
        static char const * short_name();
    private:
        fibonacci fibonaccci_;
    };
}

namespace cutrialdive {

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    lucas::residue_type lucas::value_mod_2(index_type n)
    {
        return fibonaccci_.value_mod_2(n);

    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    lucas::residue_type lucas::next_value_mod_2(residue_type r, index_type n)
    {
        return fibonaccci_.next_value_mod_2(r, n);
    }


    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    lucas::state_type   lucas::make_state(index_type n, residue_type d, mu_both_t mu)
    {
        auto st = fibonaccci_.make_state(n, d, mu);
        // st.previous is F_{n-1}
        // st.current is F_n
        // We use L_{n-1} = 2F_n - F_{n-1} and L_n = F_n + 2F_{n-1} to get what we want...
        auto fnm1 = st.previous;
        auto fn = st.current;
        st.previous = submod(doublemod(fn, d), fnm1, d);
        st.current = addmod(fn, doublemod(fnm1, d), d);
        return st;
    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    lucas::residue_type lucas::value_mod(index_type n, residue_type d, state_type & st)
    {
        return fibonaccci_.value_mod(n, d, st);
    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    lucas::residue_type lucas::next_value_mod(residue_type r, index_type n, residue_type d, state_type & st)
    {
        return fibonaccci_.next_value_mod(r, n, d, st);
    }
}
