/*
* Created on 2026.06.10
* Copyright (c) Youcef Lemsafer
*/
#include "proth.hpp"

#include <stdexcept>

#include "hgint.hpp"

namespace cutrialdive {

    proth_math::proth_math(uint64_t k)
        : k_(k)
    {
        if(!k) {
            throw std::invalid_argument("Proth parameter k must be non-zero");
        }
    }

    HgInt
    proth_math::value(uint64_t n)
    {
        return HgInt{k_} * pow(HgInt{2}, n) + 1;
    }

    std::ostream&
    proth_math::print_value(uint64_t n, std::ostream & out)
    {
        return out << value(n);
    }

    std::ostream&
    proth_math::print_expression(uint64_t n, std::ostream & out)
    {
        return out << k_ << "*2^" << n << " + 1";
    }

    proth::proth(proth_math const & mathSequence)
        : math_sequence_(mathSequence)
        , short_name_("Pr_" + std::to_string(math_sequence_.get_k()))
    {}

}
