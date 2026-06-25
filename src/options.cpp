/*
* MIT License
* Created on 2026.06.18
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#include "options.hpp"


namespace cutrialdive {

    runtime_options runtime_options::default_options()
    {
        return {
            .threads_count = 4,
            .verbosity_level = verbosity::normal
        };
    }


}
