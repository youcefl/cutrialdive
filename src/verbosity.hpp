/*
* MIT License
* Created on 2026.06.18
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#pragma once

#include <cstdint>
#include <type_traits>

namespace cutrialdive {

    enum class verbosity : uint16_t
    {
        silent = 0,
        normal = 1,
        verbose = 2,
        debug   = 3,
        diagnostic = 4
    };

}

namespace cutrialdive {

    inline bool is_valid_verbosity(uint16_t verbosityLevel)
    {
        // (verbosity::silent <= verbosityLevel) is a pointless comparison of an unsigned int with zero
        return (verbosityLevel <= static_cast<std::underlying_type_t<verbosity>>(verbosity::diagnostic));
    }

}
