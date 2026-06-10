/*
* Created on 2026.05.04
* Copyright (c) Youcef Lemsafer
*/
#pragma once

#include <atomic>
#include <cstdint>
#include <vector>
#include <ostream>
#include <istream>
#include "factoring_results.hpp"

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

    private:
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
            auto count = factors_count_[i];
            if (count > max_factors_per_number_) {
                // @todo: report maximum number of factors per number exceeded while processing S(i)
                count = max_factors_per_number_;
            }
            auto span = std::span{&factors_[i * max_factors_per_number_], count};
            sort_factors(span);
            results.add_factors(span);
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
}
