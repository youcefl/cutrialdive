/*
* Creation date: 2026.04.29
* Created by Youcef Lemsafer
*/
#pragma once

#include <ostream>

#include "common_defines.h"
#include "hgint.hpp"
#include "barrett_mu_types.hpp"
#include "modular_arithmetic_detail.hpp"


namespace cutrialdive {

    // Mersenne numbers, M(n) = 2^n - 1
    struct mersenne
    {
        using index_type = uint64_t;
        using value_type = HgInt;
        using residue_type = uint64_t;
        static constexpr bool initialize_from_value = false;

        static char const* short_name();
        static value_type value(index_type n);
        static std::ostream & print_value(index_type n, std::ostream & out);
        static std::ostream & print_expression(index_type n, std::ostream & out);
        CUTRIALDIVE_DEVICE_AND_HOST static residue_type value_mod_2(index_type n);
        CUTRIALDIVE_DEVICE_AND_HOST static residue_type value_mod_mu(index_type n, residue_type d, barrett_mu_both_t mu);
        CUTRIALDIVE_DEVICE_AND_HOST static residue_type next_value_mod(residue_type v_n_mod_p, index_type n, residue_type p);
        CUTRIALDIVE_DEVICE_AND_HOST static residue_type next_value_mod_2(residue_type v_n_mod_p, index_type n);
    };
}

namespace cutrialdive {

    inline
    char const*
    mersenne::short_name()
    {
        return "M";
    }

    inline
    typename mersenne::value_type
    mersenne::value(index_type n)
    {
        return pow(value_type{2}, n) - 1;
    }

    inline
    std::ostream &
    mersenne::print_value(index_type n, std::ostream & out)
    {
        return out << value(n);
    }

    inline
    std::ostream &
    mersenne::print_expression(index_type n, std::ostream & out)
    {
        return out << "2^" << n << "-1";
    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    uint64_t mersenne::value_mod_2(
        uint64_t n
    )
    {
        return n > 0;
    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    uint64_t mersenne::value_mod_mu(
        uint64_t n,
        uint64_t d,
        barrett_mu_both_t mu
    )
    {
        if(n <= 1) {
            return ((n == 0) || (d == 1)) ? 0 : 1;
        }
        residue_type r = 1;
        residue_type m = 2;
        for(; n; n >>= 1) {
            if(n & 1) {
                r = mulmod(r, m, d, mu.mu64, mu.mu128);
            }
            m = mulmod(m, m, d, mu.mu64, mu.mu128);
        }
        return r ? r - 1 : d - 1;
    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    uint64_t mersenne::next_value_mod_2(
        uint64_t v_n_mod_p,
        uint64_t n
    )
    {
        (void)v_n_mod_p, (void)n;
        // 2^n-1 is odd when n > 0 so just return one and be done with it
        return 1;
    }


    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    uint64_t mersenne::next_value_mod(
        uint64_t v_n_mod_d,
        uint64_t n,
        uint64_t d
    )
    {
        (void)n; // unused
        // M(n+1) = 2 * M(n) + 1
        if(v_n_mod_d < d - v_n_mod_d) {
            auto r  = (v_n_mod_d << 1) + 1;
            return r == d ? 0 : r;
        }
        return v_n_mod_d - (d - v_n_mod_d) + 1;
    }
}
