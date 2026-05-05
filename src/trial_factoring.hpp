/*
* Created on 2026.05.04
* Copyright (c) Youcef Lemsafer
*/
#pragma once

#include "mode.hpp"
#include "trial_factoring_options.hpp"

namespace cutrialdive {

    /// @brief Performs trial factoring according to given options
    /// @param opts trial factoring options
    void trial_factor(
        mode_flag modeFlag,
        trial_factoring_options const & opts
    );

}
