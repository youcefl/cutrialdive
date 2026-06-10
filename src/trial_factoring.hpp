/*
* Created on 2026.05.04
* Copyright (c) Youcef Lemsafer
*/
#pragma once

#include <ostream>
#include <filesystem>

#include "num_seq_spec.hpp"
#include "trial_factoring_options.hpp"
#include "factoring_results.hpp"

namespace cutrialdive {

    /// @brief Performs trial factoring according to given options
    /// @param opts trial factoring options
    factoring_results<uint64_t, uint32_t>
    trial_factor(
        num_seq_spec numSeqSpec,
        trial_factoring_options const & opts,
        tf_runtime_options const & runtimeOpts,
        std::ostream & out
    );

    /// Resumes a trial factoring
    /// @param[in] checkpointPath path of the checkpoint file
    /// @param[in/out] out stream where to write user messages
    factoring_results<uint64_t, uint32_t>
    resume_trial_factoring(
        std::filesystem::path const & checkpointPath,
        tf_runtime_options const & runtimeOpts,
        std::ostream & out
    );
}
