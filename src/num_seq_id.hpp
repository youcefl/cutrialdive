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
        smarandache
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
        }
        return "<unknown number sequence id>";
    }

}
