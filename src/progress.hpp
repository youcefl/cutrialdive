/*
* MIT License
* Created on 2026.06.05
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#pragma once

#include <cstdint>
#include <ostream>
#include <chrono>

namespace cutrialdive {

    class progress
    {
    public:
        /// Constructs an instance
        /// @pre @c start < @c target
        /// @param[in] start first value in the range to process
        /// @param[in] target target value
        /// @param[in] period minimum duration between two refreshes of the progress
        /// @param[in/out] out stream to write the progress to
        progress(uint64_t start, uint64_t target, std::chrono::milliseconds period, std::ostream & out);
        ~progress();
        uint64_t start() const;
        uint64_t target() const;
        /// @pre current <= target()
        void update(uint64_t current, uint64_t lastPrime);
        /// Signals the end to this instance
        void end();

    private:
        uint64_t start_;
        uint64_t target_;
        std::chrono::milliseconds period_;
        std::ostream & out_;
        bool ended_;
        size_t target_digits_;
        decltype(std::chrono::steady_clock::now()) start_time_;
        decltype(std::chrono::steady_clock::now()) last_update_time_;
    };

}


// Impl
namespace cutrialdive {

    inline
    uint64_t progress::start() const
    {
        return start_;
    }

    inline
    uint64_t progress::target() const
    {
        return target_;
    }

}
