/*
* Creation date: 2026.04.29
* Created by Youcef Lemsafer
*/
#pragma once

#include <ostream>

#include "hgint.hpp"
#include "number_sequence.hpp"
#include "smarandache.hpp"

namespace cutrialdive {

    template <>
    struct number_sequence<mode_flag::smarandache>
    {
        using index_type = uint64_t;
        using value_type = HgInt;
        using residue_type = uint64_t;
        static char const* short_name();
        static value_type value(index_type n);
        static void print_value(index_type n, std::ostream & out);
        static void print_expression(index_type n, std::ostream & out);
        static residue_type next_value_mod_mu(residue_type v_n_mod_p, index_type n, residue_type p, residue_type mu_p);
        static residue_type next_value_mod_2(residue_type v_n_mod_p, index_type n);
    };

    template <>
    struct number_sequence<mode_flag::mersenne>
    {
        using index_type = uint64_t;
        using value_type = HgInt;
        using residue_type = uint64_t;
        static char const* short_name();
        static value_type value(index_type n);
        static void print_value(index_type n, std::ostream & out);
        static void print_expression(index_type n, std::ostream & out);
        static residue_type next_value_mod_mu(residue_type v_n_mod_p, index_type n, residue_type p, residue_type mu_p);
        static residue_type next_value_mod_2(residue_type v_n_mod_p, index_type n);
    };
}

namespace cutrialdive {

    // Smarandache
    inline
    char const*
    number_sequence<mode_flag::smarandache>::short_name()
    {
        return "Sm";
    }

    inline
    typename number_sequence<mode_flag::smarandache>::value_type
    number_sequence<mode_flag::smarandache>::value(index_type n)
    {
        return smarandache<value_type>(n);
    }

    inline
    void
    number_sequence<mode_flag::smarandache>::print_value(index_type n, std::ostream & out)
    {
        output_smarandache(n, out);
    }

    inline
    void
    number_sequence<mode_flag::smarandache>::print_expression(index_type n, std::ostream & out)
    {
        output_smarandache_expression(n, out);
    }

    inline
    uint64_t number_sequence<mode_flag::smarandache>::next_value_mod_mu(
        uint64_t v_n_mod_p,
        uint64_t n,
        uint64_t p,
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
    
    inline
    uint64_t number_sequence<mode_flag::smarandache>::next_value_mod_2(
        uint64_t v_n_mod_p,
        uint64_t n
    )
    {
        return (n & 1) ^ 1;
    }


    // Mersenne
    inline
    char const*
    number_sequence<mode_flag::mersenne>::short_name()
    {
        return "M";
    }

    inline
    typename number_sequence<mode_flag::mersenne>::value_type
    number_sequence<mode_flag::mersenne>::value(index_type n)
    {
        return pow(value_type{2}, n) - 1;
    }

    inline
    void
    number_sequence<mode_flag::mersenne>::print_value(index_type n, std::ostream & out)
    {
        //@todo: operator<<(std::ostream&) for HgInt
        out << "not implemented yet";
    }

    inline
    void
    number_sequence<mode_flag::mersenne>::print_expression(index_type n, std::ostream & out)
    {
        out << "2^" << n << "-1";
    }

    inline
    uint64_t number_sequence<mode_flag::mersenne>::next_value_mod_mu(
        uint64_t v_n_mod_p,
        uint64_t n,
        uint64_t p,
        uint64_t mu_p
    )
    {
        (void)n; // unused
        // M(n+1) = 2 * M(n) + 1
        auto a = (__uint128_t(v_n_mod_p) << 1) + 1;
        auto q = (a * mu_p) >> 64;
        auto r =  a - q * p;
        return r >= p ? r - p : r;
    }

    inline
    uint64_t number_sequence<mode_flag::mersenne>::next_value_mod_2(
        uint64_t v_n_mod_p,
        uint64_t n
    )
    {
        (void)v_n_mod_p, (void)n;
        // 2^n-1 is odd when n > 0 so just return one and be done with it
        return 1;
    }

}
