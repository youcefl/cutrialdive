/*
* Creation date: 2026.04.29
* Created by Youcef Lemsafer
*/
#pragma once

#include <string>

namespace cutrialdive {

    /// Available modes
    enum class mode_flag {
        mersenne,
        smarandache
    }

    /// @brief Returns the mode as a cutrialdive::mode value
    mode_flag mode_from_string(std::string const & modeAsStr);

    template <mode_flag flag>
    struct mode {
        using number_sequence_type = number_sequence<flag>;
    };
}

