/*
* MIT License
* Created on 2026.05.03
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#pragma once

#include <cstdint>

#include "common_defines.h"
#include "device_span.hpp"
#include "device_vector.hpp"
#include "factors_buffer.hpp"
#include "factoring_results.hpp"


namespace cutrialdive {

    /// Holds the factors found by the device
    template <typename ValueT>
    class device_factors_buffer
    {
    public:
        device_factors_buffer(uint64_t n0, uint64_t numbersCount, uint32_t maxFactorsPerNumber);
        /// Constructs an instance from a host factors_buffer
        template <typename PrimeT>
        device_factors_buffer(factors_buffer<PrimeT> const & factorsBuf);
        
        struct device_view_t {
            device_span<uint32_t> counts;
            device_span<ValueT> factors;
            uint32_t const max_factors;

            __device__ void push_factor(uint64_t index, ValueT factor);
        };

        device_view_t device_view() const;

        template <typename PrimeT>
        void to_host_factors_buffer(factors_buffer<PrimeT> & hostFactorsBuf) const;
        /// @overload
        template <typename PrimeT>
        factors_buffer<PrimeT> to_host_factors_buffer() const;

    private:
        device_vector<uint32_t> factors_count_;
        device_vector<uint64_t> factors_;
        uint32_t const max_factors_per_number_;
        uint64_t numbers_count_;
        uint64_t n0_;
    };
}

namespace cutrialdive {

    template <typename ValueT>
    inline
    device_factors_buffer<ValueT>::device_factors_buffer(uint64_t n0, uint64_t numbersCount, uint32_t maxFactorsPerNumber)
        : max_factors_per_number_(maxFactorsPerNumber)
        , numbers_count_(numbersCount)
        , n0_(n0)
    {
        factors_count_.assign_zero(numbersCount);
        factors_.resize(max_factors_per_number_ * numbersCount);
    }

    template <typename ValueT>
    template <typename PrimeT>
    inline
    device_factors_buffer<ValueT>::device_factors_buffer<PrimeT>(
        factors_buffer<PrimeT> const & factorsBuf
      ) : max_factors_per_number_(factorsBuf.max_factors_per_number_)
        , numbers_count_(factorsBuf.numbers_count_)
        , n0_(factorsBuf.n0_)
    {
        factors_count_.resize(factorsBuf.factors_count_.size());
        copy_to_device(factors_count_.data(), &factorsBuf.factors_count_[0], factors_count_.size());
        factors_.resize(factorsBuf.factors_.size());
        copy_to_device(factors_.data(), & factorsBuf.factors_[0], factors_.size());
    }

    template <typename ValueT>
    template <typename PrimeT>
    inline
    factors_buffer<PrimeT> device_factors_buffer<ValueT>::to_host_factors_buffer() const
    {
        static_assert(std::is_same_v<PrimeT, ValueT>); // PrimeT != ValueT not supported
        factors_buffer<PrimeT> factorsBuf{n0_, numbers_count_, max_factors_per_number_};
        copy_from_device(&factorsBuf.factors_count_[0], factors_count_.data(), factors_count_.size());
        copy_from_device(&factorsBuf.factors_[0], factors_.data(), factors_.size());
        return std::move(factorsBuf);
    }

    template <typename ValueT>
    template <typename PrimeT>
    inline
    void device_factors_buffer<ValueT>::to_host_factors_buffer(factors_buffer<PrimeT> & hostFactorsBuf) const
    {
        static_assert(std::is_same_v<PrimeT, ValueT>); // PrimeT != ValueT not supported

        if( (n0_ == hostFactorsBuf.n0_)
          && (numbers_count_ == hostFactorsBuf.numbers_count_)
          && (max_factors_per_number_ == hostFactorsBuf.max_factors_per_number_) ) {
            copy_from_device(&hostFactorsBuf.factors_count_[0], factors_count_.data(), factors_count_.size());
            copy_from_device(&hostFactorsBuf.factors_[0], factors_.data(), factors_.size());
            return;
        }
        hostFactorsBuf = to_host_factors_buffer<PrimeT>();
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

}
