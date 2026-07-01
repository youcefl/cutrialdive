/*
* MIT License
* Created on 2026.06.10
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#pragma once

#include <cstdint>
#include <ostream>
#include <stdexcept>

#include "plus_or_minus.hpp"
#include "common_defines.h"
#include "hgint.hpp"
#include "barrett_mu_types.hpp"
#include "modular_arithmetic_detail.hpp"

namespace cutrialdive {

    /// Proth: k*2^n + 1, k > 0
    /// Riesel: k*2^n - 1, k > 0
    template <plus_or_minus PlusOrMinus>
    class riesel_proth_math
    {
    public:
        using index_type = uint64_t;
        using residue_type = uint64_t;

        /// @pre k > 0
        riesel_proth_math(uint64_t k);

        uint64_t get_k() const;

        HgInt value(uint64_t n);
        std::ostream& print_value(uint64_t n, std::ostream & out);
        std::ostream& print_expression(uint64_t n, std::ostream & out);
        CUTRIALDIVE_DEVICE_AND_HOST uint64_t value_mod_2(uint64_t n);
        CUTRIALDIVE_DEVICE_AND_HOST uint64_t value_mod_mu(uint64_t n, uint64_t d, mu_both_t mu);
        CUTRIALDIVE_DEVICE_AND_HOST uint64_t next_value_mod_2(uint64_t r, uint64_t n);
        CUTRIALDIVE_DEVICE_AND_HOST uint64_t next_value_mod(uint64_t r, uint64_t n, uint64_t d);

    private:
        uint64_t k_;
    };

    template <plus_or_minus PlusOrMinus>
    class riesel_proth
    {
    public:
        riesel_proth(riesel_proth_math<PlusOrMinus> const & mathSequence);

        riesel_proth_math<PlusOrMinus> math_sequence() const;

        char const* short_name() const;

    private:
        riesel_proth_math<PlusOrMinus> math_sequence_;
        std::string short_name_;
    };
}


namespace cutrialdive {

    template <plus_or_minus PlusOrMinus>
    inline
    uint64_t riesel_proth_math<PlusOrMinus>::get_k() const
    {
        return k_;
    }

    template <plus_or_minus PlusOrMinus>
    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    uint64_t riesel_proth_math<PlusOrMinus>::value_mod_2(
        uint64_t n
    )
    {
        return (n > 0) ? 1 : (k_ & 1) ^ 1;
    }

    template <plus_or_minus PlusOrMinus>
    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    uint64_t riesel_proth_math<PlusOrMinus>::next_value_mod_2(uint64_t, uint64_t)
    {
        return 1;
    }

    template <plus_or_minus PlusOrMinus>
    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    uint64_t riesel_proth_math<PlusOrMinus>::value_mod_mu(uint64_t n, uint64_t d, mu_both_t mu)
    {
        auto r = modpow(2, n, d, mu.mu64, mu.mu128);
        r = mulmod(r, k_, d, mu.mu64, mu.mu128);
        if constexpr(PlusOrMinus == plus_or_minus::plus) {
            return (r == d - 1) ? 0 : r + 1;
        } else if constexpr(PlusOrMinus == plus_or_minus::minus) {
            return r ? r - 1 : d - 1;
        }
    }

    template <plus_or_minus PlusOrMinus>
    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    uint64_t riesel_proth_math<PlusOrMinus>::next_value_mod(uint64_t r, uint64_t /*n*/, uint64_t d)
    {
        // proth: r = (2r - 1) mod d
        // riesel: r = (2r + 1) mod d
        r = (r < d - r) ? r << 1 : r - (d - r);
        if constexpr(PlusOrMinus == plus_or_minus::plus) {
            return r ? r - 1 : d - 1;
        } else if constexpr(PlusOrMinus == plus_or_minus::minus) {
            return (r == d - 1) ? 0 : r + 1;
        }
    }

    template <plus_or_minus PlusOrMinus>
    inline riesel_proth_math<PlusOrMinus>
    riesel_proth<PlusOrMinus>::math_sequence() const
    {
        return math_sequence_;
    }

    template <plus_or_minus PlusOrMinus>
    inline char const *
    riesel_proth<PlusOrMinus>::short_name() const
    {
        return short_name_.data();
    }

    template <plus_or_minus PlusOrMinus>
    inline
    riesel_proth_math<PlusOrMinus>::riesel_proth_math(uint64_t k)
        : k_(k)
    {
        if(!k) {
            static const auto msg = std::string{(PlusOrMinus == plus_or_minus::plus) ? "Proth" : "Riesel"}
                                        + " parameter k must be non-zero";
            throw std::invalid_argument(msg);
        }
    }

    template <plus_or_minus PlusOrMinus>
    inline
    HgInt
    riesel_proth_math<PlusOrMinus>::value(uint64_t n)
    {
        auto v = HgInt{k_} * pow(HgInt{2}, n);
        return  (PlusOrMinus == plus_or_minus::plus) ? v + 1 : v - 1;
    }

    template <plus_or_minus PlusOrMinus>
    inline
    std::ostream&
    riesel_proth_math<PlusOrMinus>::print_value(uint64_t n, std::ostream & out)
    {
        return out << value(n);
    }

    template <plus_or_minus PlusOrMinus>
    inline
    std::ostream&
    riesel_proth_math<PlusOrMinus>::print_expression(uint64_t n, std::ostream & out)
    {
        return out << k_ << "*2^" << n << " " << ((PlusOrMinus == plus_or_minus::plus) ? "+" : "-") << " 1";
    }

    template <plus_or_minus PlusOrMinus>
    inline
    riesel_proth<PlusOrMinus>::riesel_proth(riesel_proth_math<PlusOrMinus> const & mathSequence)
        : math_sequence_(mathSequence)
        , short_name_(((PlusOrMinus == plus_or_minus::plus) ? "Pr_" : "R_") + std::to_string(math_sequence_.get_k()))
    {}


    using proth = riesel_proth<plus_or_minus::plus>;
    using riesel = riesel_proth<plus_or_minus::minus>;

}
