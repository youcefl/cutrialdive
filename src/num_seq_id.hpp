/*
* MIT License
* Created on 2026.05.05
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#pragma once

#include <string>
#include <cstdint>
#include <stdexcept>

namespace cutrialdive {

    /// Available number sequence ids
    enum class num_seq_id : uint32_t
    {
        mersenne,
        proth,
        riesel,
        smarandache,
        fibonacci,
        lucas,
        proth_fixed_n,
        riesel_fixed_n,
        tribonacci
    };

    /// Returns a string
    template <typename StrT>
    StrT to_string(num_seq_id numSeqId);

    /// @brief Returns the number sequence id as a num_seq_id value
    num_seq_id num_seq_id_from_string(std::string const & numSeqIdAsStr);
}

namespace cutrialdive {
    namespace details {
        std::runtime_error unexpected_num_seq_id_exception(num_seq_id numSeqId);
    }

    template <typename StrT>
    StrT to_string(num_seq_id numSeqId)
    {
        switch(numSeqId) {
            case num_seq_id::mersenne: return {"mersenne"};
            case num_seq_id::proth: return {"proth"};
            case num_seq_id::riesel: return {"riesel"};
            case num_seq_id::smarandache: return {"smarandache"};
            case num_seq_id::fibonacci: return {"fibonacci"};
            case num_seq_id::lucas: return {"lucas"};
            case num_seq_id::proth_fixed_n: return {"proth_fixed_n"};
            case num_seq_id::riesel_fixed_n: return {"riesel_fixed_n"};
            case num_seq_id::tribonacci: return {"tribonacci"};
        }
        return "<unknown number sequence id>";
    }

}
