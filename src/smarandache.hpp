/*
 * Creation date: 2026.04.21
 * Created by Youcef Lemsafer
 */
#pragma once

#include <cstdint>
#include <iostream>
#include <vector>
#include <functional>
#include <optional>
#include <filesystem>
#include <chrono>


#include "factoring_results.hpp"

namespace cutrialdive
{
    /// @brief Returns Sm(n) = 1234567...n
    /// @tparam NumberT number type
    /// @param n 
    /// @return 
    template <typename NumberT>
    NumberT smarandache(uint64_t n);

    /// @brief Outputs Sm(n)
    /// @param[in/out] out the output stream
    /// @param[in] n indice
    /// @return output stream
    std::ostream &output_smarandache(std::ostream &out, uint64_t n);

    /// @brief Outputs Sm(n) as an expression
    /// @param[in/out] out the output stream
    /// @param[in] n indice
    /// @return output stream
    std::ostream &output_smarandache_expression(std::ostream &out, uint64_t n);

    struct trial_factoring_options
    {
        /// @brief Lower bound of the set of indices
        uint64_t n0;
        /// @brief Upper bound of the set of indices
        uint64_t n1;
        /// @brief Lower bound on the set of primes
        uint64_t f0;
        /// @brief Upper bound on the set of primes
        uint64_t f1;
        /// @brief Optional: path of the output file
        std::optional<std::filesystem::path> output_path;
    };

    /// @brief Performs trial factoring according to given options
    /// @param opts trial factoring options
    void trial_factor(trial_factoring_options const & opts);

    /// @brief Performs trial factoring of Sm(k) for all k in [n0, n1[
    /// using the prime numbers in [f0, f1[
    /// @param results[in/out] previous factorization results if any (pass empty vector on first run)
    void trial_factor(uint64_t n0, uint64_t n1,
        uint64_t f0, uint64_t f1,
        factoring_results<uint64_t, uint32_t> & results
    );

    /// @brief Computes the cofactor i.e. number / (product of factors)
    /// @tparam NumT number type
    /// @tparam Factor type
    /// @param number
    /// @param factors 
    template <typename NumT, typename FactorT, typename ExponentT>
    NumT compute_cofactor(NumT number, std::vector<factor<FactorT, ExponentT>> const & factors);

    /// @brief Run a PRP test on (Sm(n) / F) where F is the product of the provided factors
    /// If @param haveToBoostFactors is true we divide by each prime factor p as many times as possible
    /// and the corresponding exponent is updated, otherwise we use the factors as provided.
    /// @param n 
    /// @param factors known prime factors of Sm(n)
    /// @param haveToBoostFactors whether to boost the factors i.e. compute largest k that p^k divides Sm(n)
    void run_prp_test(uint64_t n, std::vector<factor<uint64_t, uint32_t>> & factors, bool haveToBoostFactors);
}


namespace cutrialdive {

    namespace details {
        template <typename T>
        inline auto floor_log10(T n)
        {
            auto log10n = 0;
            for(T d{1}; (d *= 10) <= n;) {
                ++log10n;
            }
            return log10n;
        }
    }

    template <typename T>
    inline T smarandache(uint64_t n)
    {
        if(n <= 1) {
            return T{n};
        }
        auto log10n = details::floor_log10(n);
        static const auto c3 = T{1490840987654358ull} * pow(T{10}, 180) - T{10990220};
        static const auto c4 = pow(T{10}, 2701) * c3 * T{10201} - T{1109890222000ull};
        static const auto c5 = pow(T{10}, 36000) * c4 * T{123454321} - T{123443211358} * pow(T{33300}, 2);
        static const std::array<std::function<T(uint64_t)>, 5> sm = {
              [](auto x) { return (pow(T{10}, x + 1) - T{9 * x + 10}) / T{81}; }
            , [](auto x) { return (pow(T{10}, 2*x - 18) * T{uint64_t(123456789)*99*99 + 991} 
                                                    - T{99 * x + 100}) / T{9801}; }
            , [&](auto x) { return (c3 * pow(T{10}, 3*x - 296) - T{120879* x + 121000}) 
                                            / T{120758121}; }
            , [&](auto x) { return (c4 * pow(T{10}, 4*(x - 999)) - T{12321} * T{ 9999*x + 10000 }) 
                                            / T{1231853592321ull}; }
            , [&](auto x) { return (c5 * pow(T{10}, 5*(x - 9999)) - T{12321} * T{1234321} * (T{99999} * T{x} + T{100000})) 
                                            / (T{12321} * T{1234321} * T{99999} * T{99999}); }
            };
        if(log10n < sm.size()) {
            return sm[log10n](n);
        } else {
            throw std::runtime_error{"smarandache(n) for n >= 10^5 not supported yet"};
        }
    }

    template <typename Func>
    inline void time(char const* prefix, Func f)
    {
        auto start_time = std::chrono::high_resolution_clock::now();
        f();
        auto end_time = std::chrono::high_resolution_clock::now();
        std::cout << prefix << std::chrono::duration<double>(end_time - start_time) << std::endl;
    }

    /// @brief Computes the cofactor by dividing @param number by the product of provided @param factors
    /// @tparam NumT Number type
    /// @tparam FactorT Factor type
    /// @param number the number whose cofactor is to be computed
    /// @param factors list of factors (p_i, e_i)
    /// @return number / (product of p_i ^ e_i)
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

    /// @brief Computes the factor by dividing @param number by the product of provided @param factors
    /// while computing the maximal exponent for each factor
    /// @tparam NumT Number type
    /// @tparam FactorT Factor type
    /// @param number the number whose cofactor is to be computed
    /// @param factors list of factors (p_i, e_i) where the {e_i} are not necessarily maximal
    /// @return (number / (product of p_i ^ k_i)) where k_i >= e_i is such that the cofactor is not divisible by p_i
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
