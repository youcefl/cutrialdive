/*
* Creation date: 2026.04.29
* Created by Youcef Lemsafer
*/
#pragma once

#include <ostream>

#include "common_defines.h"
#include "hgint.hpp"
#include "number_sequence.hpp"
#include "barrett_mu_types.hpp"
#include "modular_arithmetic_detail.hpp"


namespace cutrialdive {

    // Mersenne numbers: 2^n-1, n integer > 0
    template <>
    struct number_sequence<num_seq_id::mersenne>
    {
        using index_type = uint64_t;
        using value_type = HgInt;
        using residue_type = uint64_t;
        using barrett_mu_type = no_barrett_t;

        static char const* short_name();
        static value_type value(index_type n);
        static void print_value(index_type n, std::ostream & out);
        static void print_expression(index_type n, std::ostream & out);
        CUTRIALDIVE_DEVICE_AND_HOST static residue_type value_mod_2(index_type n);
        CUTRIALDIVE_DEVICE_AND_HOST static residue_type next_value_mod(residue_type v_n_mod_p, index_type n, residue_type p);
        CUTRIALDIVE_DEVICE_AND_HOST static residue_type next_value_mod_mu(residue_type v_n_mod_p, index_type n, residue_type p, barrett_mu_type mu_p);
        CUTRIALDIVE_DEVICE_AND_HOST static residue_type next_value_mod_2(residue_type v_n_mod_p, index_type n);
    };
}

namespace cutrialdive {

    inline
    char const*
    number_sequence<num_seq_id::mersenne>::short_name()
    {
        return "M";
    }

    inline
    typename number_sequence<num_seq_id::mersenne>::value_type
    number_sequence<num_seq_id::mersenne>::value(index_type n)
    {
        return pow(value_type{2}, n) - 1;
    }

    inline
    void
    number_sequence<num_seq_id::mersenne>::print_value(index_type n, std::ostream & out)
    {
        //@todo: operator<<(std::ostream&) for HgInt
        out << "not implemented yet";
    }

    inline
    void
    number_sequence<num_seq_id::mersenne>::print_expression(index_type n, std::ostream & out)
    {
        out << "2^" << n << "-1";
    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    uint64_t number_sequence<num_seq_id::mersenne>::value_mod_2(
        uint64_t n
    )
    {
        return n > 0;
    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    uint64_t number_sequence<num_seq_id::mersenne>::next_value_mod_2(
        uint64_t v_n_mod_p,
        uint64_t n
    )
    {
        (void)v_n_mod_p, (void)n;
        // 2^n-1 is odd when n > 0 so just return one and be done with it
        return 1;
    }


    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    uint64_t number_sequence<num_seq_id::mersenne>::next_value_mod(
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
