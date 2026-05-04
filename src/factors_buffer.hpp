/*
* Created on 2026.05.03
* Copyright (c) Youcef Lemsafer
*/
#pragma once

#include <atomic>
#include <cstdint>
#include <vector>
#include "factoring_results.hpp"

namespace cutrialdive {

    template <typename PrimeT, typename ExponentT>
    class factors_buffer {
    public:
        factors_buffer(uint64_t n0, uint64_t numbersCount, uint32_t maxFactorsPerNumber);
        factors_buffer(factors_buffer&&) = default;
        factors_buffer(factors_buffer const &) = delete;
        factors_buffer& operator=(factors_buffer const &) = delete;

        /// @brief Pushes a new factor for number of index numberIdx
        /// @pre numberIdx < numbersCount
        void push_factor(uint64_t numberIdx, PrimeT prime);

        factoring_results<PrimeT, ExponentT> to_factoring_results();

    private:
        std::unique_ptr<std::atomic_uint32_t[]> factors_count_;
        std::vector<factor<PrimeT, ExponentT>> factors_;
        uint32_t max_factors_per_number_;
        uint64_t numbers_count_;
        uint64_t n0_;
    };

    template <typename PrimeT, typename ExponentT>
    inline
    factors_buffer<PrimeT, ExponentT>::factors_buffer(uint64_t n0, uint64_t numbersCount, uint32_t maxFactorsPerNumber)
        : max_factors_per_number_(maxFactorsPerNumber)
        , numbers_count_(numbersCount)
        , n0_(n0)
    {
        factors_count_ = std::make_unique<std::atomic_uint32_t[]>(numbersCount);
        for(size_t i = 0; i < numbers_count_; ++i) {
            factors_count_[i] = uint32_t{};
        }
        factors_.assign(numbers_count_ * max_factors_per_number_, factor<PrimeT, ExponentT>{});
    }

    template <typename PrimeT, typename ExponentT>
    inline
    void factors_buffer<PrimeT, ExponentT>::push_factor(uint64_t numberIdx, PrimeT prime)
    {
        auto pos = factors_count_[numberIdx].fetch_add(1); // atomic increment!
        if(pos < max_factors_per_number_) {
            factors_[numberIdx * max_factors_per_number_ + pos] = factor<PrimeT, ExponentT>{prime, ExponentT{1}};
        }
        // we do not report inability to add the factor here, it is done later
        // by comparing factors_count_[numberIdx] to max_factors_per_number_
    }

    template <typename PrimeT, typename ExponentT>
    inline
    factoring_results<PrimeT, ExponentT> factors_buffer<PrimeT, ExponentT>::to_factoring_results()
    {
        factoring_results<PrimeT, ExponentT> results{n0_, numbers_count_};
        for(std::size_t i = 0; i < numbers_count_; ++i) {
            auto count = factors_count_[i].load();
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

}
