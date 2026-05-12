/*
* Created on 2026.05.04
* Copyright (c) Youcef Lemsafer
*/
#pragma once

#include <cstdint>
#include <ostream>
#include <optional>
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
    };

    std::ostream & operator<<(std::ostream & out, trial_factoring_options const & options);
}
