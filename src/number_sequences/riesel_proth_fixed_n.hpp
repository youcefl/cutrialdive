/*
* MIT License
* Created on 2026.06.10
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#pragma once

#include <cstdint>
#include <stdexcept>
#include "common_defines.h"
#include "plus_or_minus.hpp"
#include "hgint.hpp"


namespace cutrialdive {

    /// Proth: k*2^n + 1, k > 0
    /// Riesel: k*2^n - 1, k > 0
    template <plus_or_minus PlusOrMinus>
    class riesel_proth_fixed_n_math
    {
    public:
        using index_type = uint64_t;
        using residue_type = uint64_t;
        using mu_type = mu_both_t;
        static const bool has_state = true;
        using state_type = uint64_t;

        riesel_proth_fixed_n_math(uint64_t n);

        uint64_t get_n() const;

        /// Returns 2^n mod d
        CUTRIALDIVE_DEVICE_AND_HOST state_type make_state(index_type k, residue_type d, mu_both_t mu);

        HgInt value(index_type k);
        std::ostream& print_value(index_type k, std::ostream & out);
        std::ostream& print_expression(index_type k, std::ostream & out);
        CUTRIALDIVE_DEVICE_AND_HOST uint64_t value_mod_2(index_type k);
        CUTRIALDIVE_DEVICE_AND_HOST uint64_t value_mod_mu(index_type k, uint64_t d, mu_both_t mu, state_type & st);
        CUTRIALDIVE_DEVICE_AND_HOST uint64_t next_value_mod_2(uint64_t r, index_type k);
        CUTRIALDIVE_DEVICE_AND_HOST uint64_t next_value_mod(uint64_t r, index_type k, uint64_t d, state_type & st);

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

    using riesel_fixed_n = riesel_proth_fixed_n<plus_or_minus::minus>;
    using proth_fixed_n = riesel_proth_fixed_n<plus_or_minus::plus>;
}


namespace cutrialdive {

    template <plus_or_minus PlusOrMinus>
    inline
    riesel_proth_fixed_n_math<PlusOrMinus>::riesel_proth_fixed_n_math(uint64_t n)
        : n_(n)
    {
        if(!n) {
            throw std::invalid_argument{"Invalid n value zero for Riesel/Proth fixed n"};
        }
    }

    template <plus_or_minus PlusOrMinus>
    inline
    uint64_t riesel_proth_fixed_n_math<PlusOrMinus>::get_n() const
    {
        return n_;
    }

    template <plus_or_minus PlusOrMinus>
    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    riesel_proth_fixed_n_math<PlusOrMinus>::state_type
    riesel_proth_fixed_n_math<PlusOrMinus>::make_state(index_type k, residue_type d, mu_both_t mu)
    {
        return modpow(2, n_, d, mu.mu64, mu.mu128);
    }

    template <plus_or_minus PlusOrMinus>
    inline
    HgInt riesel_proth_fixed_n_math<PlusOrMinus>::value(index_type k)
    {
        return (PlusOrMinus == plus_or_minus::plus) ? (HgInt{k} << n_) + 1 : (HgInt{k} << n_) - 1;
    }

    template <plus_or_minus PlusOrMinus>
    inline
    std::ostream& riesel_proth_fixed_n_math<PlusOrMinus>::print_value(index_type k, std::ostream & out)
    {
        return out << value(k);
    }

    template <plus_or_minus PlusOrMinus>
    inline
    std::ostream& riesel_proth_fixed_n_math<PlusOrMinus>::print_expression(index_type k, std::ostream & out)
    {
        return out << k << "*2^" << n_ << (PlusOrMinus == plus_or_minus::plus ? '+' : '-') << "1";
    }

    template <plus_or_minus PlusOrMinus>
    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    uint64_t riesel_proth_fixed_n_math<PlusOrMinus>::value_mod_2(index_type k)
    {
        return 1;
    }

    template <plus_or_minus PlusOrMinus>
    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    uint64_t riesel_proth_fixed_n_math<PlusOrMinus>::value_mod_mu(index_type k, uint64_t d, mu_both_t mu, state_type & st)
    {
        uint64_t r = mod(k, d, mu.mu64);
        r = mulmod(r, st, d, mu.mu64, mu.mu128);
        return (PlusOrMinus == plus_or_minus::plus)
                    ? ((r == d - 1) ? 0 : r + 1)
                    : (!r ? d - 1 : r - 1);
    }

    template <plus_or_minus PlusOrMinus>
    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    uint64_t riesel_proth_fixed_n_math<PlusOrMinus>::next_value_mod_2(uint64_t r, index_type k)
    {
        return 1;
    }

    template <plus_or_minus PlusOrMinus>
    CUTRIALDIVE_DEVICE_AND_HOST uint64_t riesel_proth_fixed_n_math<PlusOrMinus>::next_value_mod(
        uint64_t r,
        index_type k,
        uint64_t d,
        state_type & st
      )
    {
        return (r >= d - st) ? r - (d - st) : r + st;
    }


    template <plus_or_minus PlusOrMinus>
    inline
    riesel_proth_fixed_n<PlusOrMinus>::riesel_proth_fixed_n(riesel_proth_fixed_n_math<PlusOrMinus> const & mathSequence)
        : math_sequence_(mathSequence)
    {
        if constexpr(PlusOrMinus == plus_or_minus::plus) {
            short_name_ = "Prn_";
        } else {
            short_name_ = "Rn_";
        }
        short_name_ += std::to_string(math_sequence_.get_n());
    }

    template <plus_or_minus PlusOrMinus>
    inline
    riesel_proth_fixed_n_math<PlusOrMinus> riesel_proth_fixed_n<PlusOrMinus>::math_sequence() const
    {
        return math_sequence_;
    }

    template <plus_or_minus PlusOrMinus>
    inline
    char const* riesel_proth_fixed_n<PlusOrMinus>::short_name() const
    {
        return short_name_.c_str();
    }

}
