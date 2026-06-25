/*
* MIT License
* Created on 2026.06.10
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#pragma once

#include <cstdint>
#include <ostream>

#include "common_defines.h"
#include "hgint.hpp"
#include "barrett_mu_types.hpp"
#include "modular_arithmetic_detail.hpp"

namespace cutrialdive {

    /// Proth numbers: k*2^n + 1, k > 0
    class proth_math
    {
    public:
        using index_type = uint64_t;
        using residue_type = uint64_t;

        /// @pre k > 0
        proth_math(uint64_t k);

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

    class proth
    {
    public:
        proth(proth_math const & mathSequence);

        proth_math math_sequence() const;

        char const* short_name() const;

    private:
        proth_math math_sequence_;
        std::string short_name_;
    };
}


namespace cutrialdive {

    inline
    uint64_t proth_math::get_k() const
    {
        return k_;
    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    uint64_t proth_math::value_mod_2(
        uint64_t n
    )
    {
        return (n > 0) ? 1 : (k_ & 1) ^ 1;
    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    uint64_t proth_math::next_value_mod_2(uint64_t, uint64_t)
    {
        return 1;
    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    uint64_t proth_math::value_mod_mu(uint64_t n, uint64_t d, mu_both_t mu)
    {
        auto r = modpow(2, n, d, mu.mu64, mu.mu128);
        r = mulmod(r, k_, d, mu.mu64, mu.mu128);
        return (r == d - 1) ? 0 : r + 1;
    }

    CUTRIALDIVE_INLINE CUTRIALDIVE_DEVICE_AND_HOST
    uint64_t proth_math::next_value_mod(uint64_t r, uint64_t /*n*/, uint64_t d)
    {
        // r = (2r - 1) mod d
        r = (r < d - r) ? r << 1 : r - (d - r);
        return r ? r - 1 : d - 1;
    }

    inline proth_math
    proth::math_sequence() const
    {
        return math_sequence_;
    }

    inline char const *
    proth::short_name() const
    {
        return short_name_.data();
    }

}
