/*
* MIT License
* Created on 2026.05.03
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#pragma once

#include <vector>
#include <span>
#include <ostream>

#include "factor.hpp"

namespace cutrialdive {

    template <typename ValueT, typename ExponentT> class factoring_results;
    template <typename ValueT, typename ExponentT>
    std::ostream & operator<<(std::ostream & out, factoring_results<ValueT, ExponentT> const & results);


    /// Results of a trial factoring run
    template <typename ValueT, typename ExponentT>
    class factoring_results {
    public:
        /// Constructs an instance with the index of the first trial factored number.
        factoring_results(uint64_t startIndex, std::size_t numbersCount);

        /// Returns the first number's index i.e. if this instances holds the results of trial factoring
        /// S(n0), ..., S(n1 - 1), then n0 is returned.
        uint64_t n0() const;

        /// Returns the number of factor sets i.e. this instance holds the factors of S(k)
        /// for all k in [n0, n0 + size()[.
        std::size_t size() const;

        /// Returns true if and only if this->size() == 0
        bool empty() const;

        /// Reserve space for @param factorsCount factors
        void reserve(std::size_t factorsCount);

        /// Adds the factors for current number
        void add_factors_with_exp(std::span<factor<ValueT, ExponentT>> factors, uint32_t numFactorsFound);
        /// @overload
        template <typename PrimeT>
        void add_factors(std::span<PrimeT> factors, uint32_t numFactorsFound);

        /// Returns the factors of the i-th number
        /// @pre i < size()
        std::span<factor<ValueT, ExponentT> const>  operator[](std::size_t i) const;

        /// Returns the maximum number of factors found for single a number
        uint32_t max_factors_count() const;

        /// Returns total number of factors found
        uint64_t factors_count() const;

        /// Writes the results to the given output stream
        friend std::ostream & operator<< <ValueT, ExponentT>(
            std::ostream &,
            factoring_results<ValueT, ExponentT> const &
        );

    private:
        /// Location of the set of factors of a number
        struct location {
            std::size_t offset;
            std::size_t length;
        };

        /// First index into the number sequence i.e. first factors are factors of S(n0_)
        uint64_t n0_;
        /// Table of (offset, length) allowing to map a TFed number to its factors
        std::vector<location> locations_;
        /// Table of the factors of all numbers processed
        std::vector<factor<ValueT, ExponentT>> factors_;
        /// Maximum number of factors found for a single sequence term (counting eventual factors not
        /// stored because of the limit max_factors_per_number_).
        uint32_t max_factors_count_;
        /// Overall number of factors found (includes eventual factors not
        /// stored because of the limit max_factors_per_number_).
        uint64_t factors_count_;
    };
}

namespace cutrialdive {

    template <typename ValueT, typename ExponentT>
    inline
    factoring_results<ValueT, ExponentT>::factoring_results(uint64_t startIndex, std::size_t numbersCount)
        : n0_(startIndex)
        , max_factors_count_{}
        , factors_count_{}
    {
        locations_.reserve(numbersCount);
    }

    template <typename ValueT, typename ExponentT>
    inline
    void factoring_results<ValueT, ExponentT>::reserve(std::size_t factorsCount)
    {
        factors_.reserve(factorsCount);
    }

    template <typename ValueT, typename ExponentT>
    inline
    uint64_t factoring_results<ValueT, ExponentT>::n0() const
    {
        return n0_;
    }

    template <typename ValueT, typename ExponentT>
    inline
    std::size_t factoring_results<ValueT, ExponentT>::size() const
    {
        return locations_.size();
    }

    template <typename ValueT, typename ExponentT>
    inline
    bool factoring_results<ValueT, ExponentT>::empty() const
    {
        return size() == 0;
    }

    template <typename ValueT, typename ExponentT>
    inline
    std::span<factor<ValueT, ExponentT> const> factoring_results<ValueT, ExponentT>::operator[](std::size_t i) const
    {
        auto const & loc = locations_[i];
        return {factors_.data() + loc.offset, loc.length};
    }

    template <typename ValueT, typename ExponentT>
    inline
    uint32_t factoring_results<ValueT, ExponentT>::max_factors_count() const
    {
        return max_factors_count_;
    }

    template <typename ValueT, typename ExponentT>
    inline
    uint64_t factoring_results<ValueT, ExponentT>::factors_count() const
    {
        return factors_count_;
    }

    template <typename ValueT, typename ExponentT>
    inline
    void factoring_results<ValueT, ExponentT>::add_factors_with_exp(
        std::span<factor<ValueT, ExponentT>> newFactors,
        uint32_t numFactorsFound
      )
    {
        auto offset = factors_.size();
        auto length = newFactors.size();
        factors_.insert(std::end(factors_), std::begin(newFactors), std::end(newFactors));
        locations_.push_back({offset, length});
        if(numFactorsFound > max_factors_count_) {
            max_factors_count_ = numFactorsFound;
        }
        factors_count_ += numFactorsFound;
    }

    template <typename ValueT, typename ExponentT>
    template <typename PrimeT>
    inline
    void factoring_results<ValueT, ExponentT>::add_factors(std::span<PrimeT> newFactors, uint32_t numFactorsFound)
    {
        auto offset = factors_.size();
        auto length = newFactors.size();
        std::for_each(std::begin(newFactors), std::end(newFactors), [&](auto const & prime){
            factors_.emplace_back(prime, ExponentT{1});
        });
        locations_.push_back({offset, length});
        if(numFactorsFound > max_factors_count_) {
            max_factors_count_ = numFactorsFound;
        }
        factors_count_ += numFactorsFound;
    }

    template <typename ValueT, typename ExponentT>
    inline
    std::ostream & operator<<(std::ostream & out, factoring_results<ValueT, ExponentT> const & results)
    {
        auto indexOfFirstNumber = results.n0_;

        for(uint64_t i = 0, iEnd = results.size(); i < iEnd; ++i) {
            auto indexOfNumber = indexOfFirstNumber + i;
            out << indexOfNumber << " - ";
            auto factors = results[i];
            char const * sep = "";
            std::for_each(std::begin(factors), std::end(factors), [&sep, &out](auto const & factor) {
                out << sep << factor;
                sep = ", ";
            });
            out << "\n";
        }
        return out;
    }
   
}