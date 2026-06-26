/*
* MIT License
* Created on 2026.05.04
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#pragma once

#include <atomic>
#include <cstdint>
#include <vector>
#include <ostream>
#include <istream>
#include <algorithm>
#include <optional>

#include "factoring_results.hpp"
#include "logger.hpp"


namespace cutrialdive {

    template <typename PrimeT>
    class factors_buffer;
    template <typename PrimeT>
    void serialize_factors_buffer(std::ostream & out, factors_buffer<PrimeT> const &);
    template <typename PrimeT>
    factors_buffer<PrimeT> deserialize_factors_buffer(std::istream & in);
    template <typename ValueT>
    class device_factors_buffer;

    /// Stores factorization data for a contiguous range of sequence terms.
    template <typename PrimeT>
    class factors_buffer {
    public:
        factors_buffer(uint64_t n0, uint64_t numbersCount, uint32_t maxFactorsPerNumber);
        factors_buffer(factors_buffer&&) noexcept = default;
        factors_buffer& operator=(factors_buffer&&) noexcept = default;
        factors_buffer(factors_buffer const &) = delete;
        factors_buffer& operator=(factors_buffer const &) = delete;

        /// @brief Pushes a new factor for the sequence term of index numberIdx
        /// @pre numberIdx < numbers_count()
        void push_factor(uint64_t numberIdx, PrimeT prime);

        template <typename ResultPrimeT, typename ExponentT>
        factoring_results<ResultPrimeT, ExponentT> to_factoring_results() const;

        /// Index of the first sequence term
        uint64_t n0() const;

        /// Returns the number of consecutive sequence values whose factors
        /// are stored in this buffer.
        ///
        /// This buffer stores factors of
        ///     S(n0_), S(n0_ + 1), ..., S(n0_ + numbers_count() - 1).
        uint64_t numbers_count() const;

        /// Returns the maximum number of factors that can be stored for a sequence term
        uint32_t max_factors_per_number() const;

        /// Checks for excess factors, if any are found reports to the user that we need
        /// more room to store all factors.
        /// O(N) in the number of sequence terms so use with caution!
        void report_excess_factors_if_any();

    private:
        /// Returns the maximum number of factors found for a single sequence term
        /// O(N) in the number of sequence terms so use with caution!
        uint32_t compute_max_factors_count();

        /// In case has_excess_factors() is true, reports to the user that we need more room to store all factors.
        /// Calls has_excess_factors() to establish whether user report has to be done.
        /// @return true if and only if report as been done.
        void report_excess_factors();

        /// Number of factors found for each number S(n0), S(n0 + 1), ..., S(n0 + numbers_count_ - 1)
        std::vector<uint32_t> factors_count_;
        /// Factors found. The factors of S(n0 + i) are at offset (i * max_factors_per_number_)
        mutable std::vector<PrimeT> factors_; // mutable because we sort the factors in to_factoring_results()
        /// The maximum number of factors that can be stored for a single number
        uint32_t max_factors_per_number_;
        /// The number of numbers this instance stores the factors for
        uint64_t numbers_count_;
        /// Index of the first sequence term this buffer stores the factors for.
        uint64_t n0_;
        /// Whether there are excess factors i.e. (compute_max_factors_count() > max_factors_per_number()).
        bool has_excess_factors_;
        /// Returns last value returned by compute_max_factors_count(), has no value if
        /// compute_max_factors_count() has not been called.
        std::optional<uint32_t> last_computed_max_factors_count_;

        friend void serialize_factors_buffer<PrimeT>(std::ostream & out, factors_buffer<PrimeT> const &);
        friend factors_buffer<PrimeT> deserialize_factors_buffer<PrimeT>(std::istream & in);
        friend class device_factors_buffer<PrimeT>;
    };

    template <typename PrimeT>
    inline
    factors_buffer<PrimeT>::factors_buffer(uint64_t n0, uint64_t numbersCount, uint32_t maxFactorsPerNumber)
        : factors_count_(numbersCount)
        , factors_(numbersCount * maxFactorsPerNumber)
        , max_factors_per_number_(maxFactorsPerNumber)
        , numbers_count_(numbersCount)
        , n0_(n0)
        , has_excess_factors_{false}
    {
    }

    template <typename PrimeT>
    inline
    void factors_buffer<PrimeT>::push_factor(uint64_t numberIdx, PrimeT prime)
    {
        auto pos = std::atomic_ref<uint32_t>{factors_count_[numberIdx]}
                        .fetch_add(1, std::memory_order_relaxed); // atomic increment!
        if(pos < max_factors_per_number_) {
            factors_[numberIdx * max_factors_per_number_ + pos] = prime;
        }
        // we do not report inability to add the factor here, it is done later
        // by comparing factors_count_[numberIdx] to max_factors_per_number_
    }

    template <typename PrimeT>
    template <typename ResultPrimeT, typename ExponentT>
    inline
    factoring_results<ResultPrimeT, ExponentT> factors_buffer<PrimeT>::to_factoring_results() const
    {
        factoring_results<ResultPrimeT, ExponentT> results{n0_, numbers_count_};

        for(std::size_t i = 0; i < numbers_count_; ++i) {
            // Can not store more than max_factors_per_number_ factors.
            auto count = (std::min)(factors_count_[i], max_factors_per_number_);
            auto span = std::span<PrimeT>{&factors_[i * max_factors_per_number_], count};
            sort_factors(span);
            results.add_factors(span, factors_count_[i]);
        }
        return results;
    }

    template <typename PrimeT>
    inline
    uint64_t factors_buffer<PrimeT>::n0() const
    {
        return n0_;
    }

    template <typename PrimeT>
    inline
    uint64_t factors_buffer<PrimeT>::numbers_count() const
    {
        return numbers_count_;
    }

    template <typename PrimeT>
    inline
    uint32_t factors_buffer<PrimeT>::max_factors_per_number() const
    {
        return max_factors_per_number_;
    }

    template <typename PrimeT>
    inline
    uint32_t factors_buffer<PrimeT>::compute_max_factors_count()
    {
        last_computed_max_factors_count_ = std::ranges::max(factors_count_);
        if(last_computed_max_factors_count_ > max_factors_per_number()) {
            has_excess_factors_ = true;
        }
        return *last_computed_max_factors_count_;
    }

    template <typename PrimeT>
    inline
    void factors_buffer<PrimeT>::report_excess_factors_if_any()
    {
        auto prevMax = last_computed_max_factors_count_.value_or(0);
        compute_max_factors_count();
        if(!has_excess_factors_
          || (*last_computed_max_factors_count_ <= prevMax)) { // report only if value has increased
            return;
        }

        auto maxFound = *last_computed_max_factors_count_;
        auto suggestedMaxFactorsPerNumber = ((maxFound / 64) + 1) * 64 + 64;

        CTDLOG_INFO()
            << "At least one number has factors count " << maxFound << " which exceeds"
            " the current storage limit (" << max_factors_per_number_ << ").\n" <<
            "Trial factoring continues normally, but some factors may be omitted from the recorded factor list for those terms.\n"
            "Consider restarting with\n"
            "  --max-factors-per-number " << suggestedMaxFactorsPerNumber << std::endl;
    }
}
