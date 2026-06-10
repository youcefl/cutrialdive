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
        /// @brief Whether to display progress
        bool is_progress_enabled{true};

        bool operator==(trial_factoring_options const&) const = default;        
    };

    std::ostream & operator<<(std::ostream & out, trial_factoring_options const & options);

    /// Execution options
    struct tf_runtime_options
    {
        std::chrono::seconds checkpoint_period;
        std::chrono::milliseconds progress_period;

        static tf_runtime_options default_options();
    };
}
