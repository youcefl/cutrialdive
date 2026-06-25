/*
* MIT License
* Created on 2026.06.18
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#pragma once

#include <cstdint>

#include "verbosity.hpp"

namespace cutrialdive {

    /// Global options
    struct runtime_options
    {
        uint32_t threads_count;
        verbosity verbosity_level;

        static runtime_options default_options();
    };

}
