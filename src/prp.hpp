/*
* Creation date: 2026.05.05
* Created by Youcef Lemsafer
*/
#pragma once

#include <cstdint>
#include <vector>
#include <algorithm>

#include "num_seq_spec.hpp"
#include "factor.hpp"


namespace cutrialdive {

    void run_prp_test(
        num_seq_spec numSeqSpec,
        uint64_t n,
        std::vector<factor<uint64_t, uint32_t>> & factors,
        bool haveToBoostFactors,
        std::ostream & out
    );

    /// @brief Computes the cofactor by dividing @param number by the product of provided @param factors
    /// @tparam NumT Number type
    /// @tparam FactorT Factor type
    /// @param number the number whose cofactor is to be computed
    /// @param factors list of factors (p_i, e_i)
    /// @return number / (product of p_i ^ e_i)
    template <typename NumT, typename FactorT, typename ExponentT>
    NumT compute_cofactor_exact(NumT number, std::vector<factor<FactorT, ExponentT>> const & factors);

    /// @brief Computes the factor by dividing @param number by the product of provided @param factors
    /// while computing the maximal exponent for each factor
    /// @tparam NumT Number type
    /// @tparam FactorT Factor type
    /// @param number the number whose cofactor is to be computed
    /// @param factors list of factors (p_i, e_i) where the {e_i} are not necessarily maximal
    /// @return (number / (product of p_i ^ k_i)) where k_i >= e_i is such that the cofactor is not divisible by p_i
    template <typename NumT, typename FactorT>
    NumT compute_cofactor_boosted(NumT number, std::vector<factor<FactorT, uint32_t>> & factors);
}


namespace cutrialdive {

    template <typename NumT, typename FactorT, typename ExponentT>
    inline NumT compute_cofactor_exact(NumT number, std::vector<factor<FactorT, ExponentT>> const & factors)
    {
        if(factors.empty()) {
            return number;
        }
        NumT product_of_factors{1};
        std::for_each(std::begin(factors), std::end(factors), 
            [&product_of_factors](auto const & factor) {
                product_of_factors *= pow(NumT{factor.prime}, factor.exponent);
        });
        return number / product_of_factors;
    }

    template <typename NumT, typename FactorT>
    inline NumT compute_cofactor_boosted(NumT number, std::vector<factor<FactorT, uint32_t>> & factors)
    {
        if(factors.empty()) {
            return number;
        }
        NumT product_of_factors{1};
        std::for_each(std::begin(factors), std::end(factors), 
            [&product_of_factors](auto const & factor) {
                product_of_factors *= pow(NumT{factor.prime}, factor.exponent);
        });
        number /= product_of_factors;
        NumT quotient, remainder;
        std::for_each(std::begin(factors), std::end(factors),
            [&](auto & factor) {
                while(true) {
                    div_mod(number, factor.prime, quotient, remainder);
                    if(!!remainder) {
                        break;
                    }
                    number = NumT{quotient};
                    ++factor.exponent;
                }
        });
        return number;
    }

}

