/*
* Created on 2026.05.04
* Copyright (c) Youcef Lemsafer
*/
#pragma once

#include <cstdint>
#include <optional>
#include <ostream>
#include <chrono>
#include <filesystem>
#include "options.hpp"

namespace cutrialdive {

    struct trial_factoring_options
    {
        /// @brief Lower bound of the set of indices
        uint64_t n0;
        /// @brief Upper bound of the set of indices
        uint64_t n1;
        /// @brief Lower bound on the set of primes
        uint64_t f0;
        /// @brief Upper bound on the set of primes
        uint64_t f1;
        /// @brief Optional: path of the output file
        std::optional<std::filesystem::path> output_path;
        /// @brief Maximum number of factors per number (default = 128)
        uint32_t max_factors_per_number{128};

        bool operator==(trial_factoring_options const&) const = default;        
    };

    std::ostream & operator<<(std::ostream & out, trial_factoring_options const & options);

    /// Execution options
    struct tf_runtime_options
    {
        /// Common options
        runtime_options common_options;
        /// @brief Checkpointing period in seconds: minimum time between checkpoint file writes.
        std::chrono::seconds checkpoint_period;
        /// @brief Progress period in milliseconds: minimum time between progress output
        std::chrono::milliseconds progress_period;
        /// @brief Whether to display progress
        bool is_progress_enabled;
        /// @brief length of sieve segments
        uint64_t segment_length;
#ifdef CUTRIALDIVE_ENABLE_GPU
        int32_t device_id;
        std::optional<int32_t> grid_size;
        std::optional<int32_t> block_size;
#endif

        static tf_runtime_options default_options();
    };
}
