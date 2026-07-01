/*
* MIT License
* Created on 2026.06.10
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#pragma once

#include <cstdint>
#include "plus_or_minus.hpp"
#include "HgInt.hpp"


namespace cutrialdive {

    /// Proth: k*2^n + 1, k > 0
    /// Riesel: k*2^n - 1, k > 0
    template <plus_or_minus PlusOrMinus>
    class riesel_proth_fixed_n_math
    {
    public:
        using index_type = __uint128_t;
        using residue_type = uint64_t;

        riesel_proth_fixed_n_math(uint64_t n);

        uint64_t get_n() const;

        HgInt value(__uint128_t k);
        std::ostream& print_value(__uint128_t k, std::ostream & out);
        std::ostream& print_expression(__uint128_t k, std::ostream & out);
        CUTRIALDIVE_DEVICE_AND_HOST uint64_t value_mod_2(__uint128_t k);
        CUTRIALDIVE_DEVICE_AND_HOST uint64_t value_mod_mu(__uint128_t k, uint64_t d, mu_both_t mu);
        CUTRIALDIVE_DEVICE_AND_HOST uint64_t next_value_mod_2(uint64_t r, __uint128_t k);
        CUTRIALDIVE_DEVICE_AND_HOST uint64_t next_value_mod(uint64_t r, __uint128_t k, uint64_t d);

    private:
        uint64_t n_;
    };

    template <plus_or_minus PlusOrMinus>
    class riesel_proth_fixed_n
    {
    public:
        riesel_proth_fixed_n(riesel_proth_fixed_n_math<PlusOrMinus> const & mathSequence);

        riesel_proth_fixed_n_math<PlusOrMinus> math_sequence() const;

        char const* short_name() const;

    private:
        riesel_proth_fixed_n_math<PlusOrMinus> math_sequence_;
        std::string short_name_;
    };

}


namespace cutrialdive {

    template <plus_or_minus PlusOrMinus>
    inline
    riesel_proth_fixed_n_math<PlusOrMinus>::riesel_proth_fixed_n_math(uint64_t n)
        : n_(n)
    {
    }

    template <plus_or_minus PlusOrMinus>
    inline
    uint64_t riesel_proth_fixed_n_math<PlusOrMinus>::get_n() const
    {
        return n_;
    }

    template <plus_or_minus PlusOrMinus>
    HgInt riesel_proth_fixed_n_math<PlusOrMinus>::value(__uint128_t k)
    {
        return (PlusOrMinus == plus_or_minus::plus) ? (HgInt{k} << n_) + 1 : (HgInt{k} << n_) - 1;
    }

    template <plus_or_minus PlusOrMinus>
    std::ostream& riesel_proth_fixed_n_math<PlusOrMinus>::print_value(__uint128_t k, std::ostream & out)
    {
        return out << value(k);
    }

    template <plus_or_minus PlusOrMinus>
    std::ostream& riesel_proth_fixed_n_math<PlusOrMinus>::print_expression(__uint128_t k, std::ostream & out)
    {

    }

    template <plus_or_minus PlusOrMinus>
    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    uint64_t riesel_proth_fixed_n_math<PlusOrMinus>::value_mod_2(__uint128_t k)
    {
        return !n_ ? (k & 1) ^ 1 : 1;
    }

    template <plus_or_minus PlusOrMinus>
    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    uint64_t riesel_proth_fixed_n_math<PlusOrMinus>::value_mod_mu(__uint128_t k, uint64_t d, mu_both_t mu)
    {
        uint64_t r;
        if(k >> 64) {
            auto q = mul128hi(k, mu.mu128);
            auto tr = k - q*k;
            r = uint64_t(tr >= d ? tr - d : tr);
        } else {
            r = mod(k, d, mu.mu64);
        }
        r = mulmod(r, modpow(2, n_, d, mu.mu64, mu.mu128), d, mu.mu64, mu.mu128);
        return (PlusOrMinus == plus_or_minus::plus)
                    ? ((r == d - 1) ? 0 : r + 1)
                    : (!r ? d - 1 : r - 1);
    }

    template <plus_or_minus PlusOrMinus>
    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    uint64_t riesel_proth_fixed_n_math<PlusOrMinus>::next_value_mod_2(uint64_t r, __uint128_t k)
    {
        return (n_ == 0) ? k_ & 1 : 1;
    }

    template <plus_or_minus PlusOrMinus>
    CUTRIALDIVE_DEVICE_AND_HOST uint64_t riesel_proth_fixed_n_math<PlusOrMinus>::next_value_mod_mu(
        uint64_t r,
        __uint128_t k,
        uint64_t d
      )
    {
        r +    
    }

}