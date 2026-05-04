/*
* Created on 2026.05.03
* Copyright (c) Youcef Lemsafer
*/
#pragma once

#include <cstdint>

#include "device_span.hpp"
#include "device_vector.hpp"
#include "factoring_results.hpp"

namespace cutrialdive {

    /// Holds the factors found by the device
    template <typename ValueT>
    class device_factors_buffer
    {
    public:
        device_factors_buffer(uint64_t n0, uint64_t numbersCount, uint32_t maxFactorsPerNumber);
        
        struct device_view_t {
            device_span<uint32_t> counts;
            device_span<ValueT> factors;
            uint32_t const max_factors;

            __device__ void push_factor(uint64_t index, ValueT factor);
        };

        device_view_t device_view() const;

        template <typename PrimeT, typename ExponentT>
        factoring_results<PrimeT, ExponentT> to_factoring_results() const;

    private:
        device_vector<uint32_t> factors_count_;
        device_vector<uint64_t> factors_;
        uint32_t const max_factors_per_number_;
        uint64_t n0_;
    };
}

namespace cutrialdive {

    template <typename ValueT>
    inline
    device_factors_buffer<ValueT>::device_factors_buffer(uint64_t n0, uint64_t numbersCount, uint32_t maxFactorsPerNumber)
        : max_factors_per_number_(maxFactorsPerNumber)
        , n0_(n0)
    {
        factors_count_.assign_zero(numbersCount);
        factors_.resize(max_factors_per_number_ * numbersCount);
    }

    template <typename ValueT>
    inline
    typename device_factors_buffer<ValueT>::device_view_t device_factors_buffer<ValueT>::device_view() const
    {
        return {factors_count_.span(), factors_.span(), max_factors_per_number_};
    }

    template <typename ValueT>
    __device__ __forceinline__
    void device_factors_buffer<ValueT>::device_view_t::push_factor(uint64_t index, ValueT factor)
    {
        auto pos = atomicAdd(&counts[index], 1);
        if(pos < max_factors) {
            factors[index * max_factors + pos] = factor;
        }
        // We do not report inability to add the factor here, it is done later
        // by comparing factors_count_[numberIdx] to max_factors_per_number_
    }

    template <typename ValueT>
    template <typename PrimeT, typename ExponentT>
    inline
    factoring_results<PrimeT, ExponentT> device_factors_buffer<ValueT>::to_factoring_results() const
    {
        factoring_results<PrimeT, ExponentT> results{n0_, factors_count_.size()};
        auto factorsCount = factors_count_.to_host();
        auto allFactors = factors_.to_host();
        std::vector<factor<PrimeT, ExponentT>> currentFactors;
        currentFactors.reserve(max_factors_per_number_);
        std::size_t i = 0;
        for(auto count : factorsCount) {
            if (count > max_factors_per_number_) {
                // @todo: report maximum number of factors per number exceeded while processing S(i)
                count = max_factors_per_number_;
            }
            currentFactors.clear();
            for(auto k = i * max_factors_per_number_; k < i * max_factors_per_number_ + count; ++k) {
                currentFactors.push_back(factor<PrimeT, ExponentT>{allFactors[k], ExponentT{1}});
            }
            auto span = std::span{std::begin(currentFactors), std::end(currentFactors)};
            sort_factors(span);
            results.add_factors(span);
            ++i;
        }
        return results;
    }

}
