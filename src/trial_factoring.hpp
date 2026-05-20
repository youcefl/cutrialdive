/*
* Created on 2026.05.04
* Copyright (c) Youcef Lemsafer
*/
#pragma once

#include "num_seq_spec.hpp"
#include "trial_factoring_options.hpp"

namespace cutrialdive {

    /// @brief Performs trial factoring according to given options
    /// @param opts trial factoring options
    void trial_factor(
        num_seq_spec numSeqSpec,
        trial_factoring_options const & opts
    );

}
