/*
* MIT License
* Created on 2026.05.19
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#pragma once

#include <cstdint>
#include <functional>
#include <array>
#include <mutex>
#include <atomic>

namespace cutrialdive {

    /// A lazily initialized thread safe array
    template <typename T, size_t N>
    class lazy_init_array
    {
    public:
        /// Constructs an instance with the value at index 0 and a function to be called
        /// when initialization i-th array cell is requested
        lazy_init_array(T&& firstValue, std::function<T(T const &, std::size_t)>&& computeNext);

        /// Returns element at index i
        /// Elements are added if needed i.e. i-th has not been initialized yet
        /// @pre i < N
        T operator[](std::size_t i) const;

    private:
        void ensure(std::size_t i) const;

        mutable std::array<T, N> array_;
        std::function<T(T const &, std::size_t)> compute_next_;
        mutable std::mutex mutex_;
        mutable std::atomic<size_t> initialized_count_{};

    };

}

namespace cutrialdive {

    template <typename T, size_t N>
    inline
    lazy_init_array<T, N>::lazy_init_array(T&& firstValue, std::function<T(T const &, std::size_t)>&& computeNext)
        : compute_next_(std::move(computeNext))
    {
        array_[0] = std::move(firstValue);
        initialized_count_.store(1, std::memory_order_release);
    }

    template <typename T, size_t N>
    inline
    T lazy_init_array<T, N>::operator[](std::size_t i) const
    {
        ensure(i);
        return array_[i];
    }

    template <typename T, size_t N>
    inline
    void lazy_init_array<T, N>::ensure(std::size_t i) const
    {
        if(i < initialized_count_.load(std::memory_order_acquire)) {
            return;
        }
        std::lock_guard lock{mutex_};
        if(i < initialized_count_.load(std::memory_order_relaxed)) {
            return;
        }
        auto last = initialized_count_.load(std::memory_order_relaxed);
        // last is always > 0 because the constructors sets initialized_count_ to 1
        for(auto j = last; j <= i; ++j) {
            array_[j] = compute_next_(array_[j - 1], j);
        }
        initialized_count_.store(i + 1, std::memory_order_release);
    }

}
