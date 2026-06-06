/*
* Created on 2026.06.05
* Copyright (c) Youcef Lemsafer
*/
#pragma once

#include <cstdint>
#include <ostream>
#include <chrono>

namespace cutrialdive {

    class progress
    {
    public:
        progress(uint64_t target, std::ostream & out);
        ~progress();
        uint64_t target() const;
        /// @pre current <= target()
        void update(uint64_t current, uint64_t lastPrime);
        void end();
    private:
        uint64_t target_;
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
    uint64_t progress::target() const
    {
        return target_;
    }
}