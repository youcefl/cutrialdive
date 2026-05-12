/*
* Creation date: 2026.05.04
* Copyright (c) Youcef Lemsafer
*/
#pragma once

#include <stdexcept>
#include <vector>
#include <string>
#include <string_view>
#include <sstream>
#include <ostream>
#include <regex>
#include <ranges>


namespace cutrialdive {

    /// A prime factor with its exponent
    template <typename PrimeT, typename ExponentT>
    struct factor {
        PrimeT prime;
        ExponentT exponent;
    };

    /// Sorts the factors in the given span (from lowest to highest)
    template <typename PrimeT, typename ExponentT>
    void sort_factors(std::span<factor<PrimeT, ExponentT>> span);

    /// Outputs a factor to the given stream
    template <typename PrimeT, typename ExponentT>
    std::ostream & operator<<(std::ostream & out, factor<PrimeT, ExponentT> const & factor);

    template <typename PrimeT, typename ExponentT>
    std::vector<factor<PrimeT, ExponentT>> parse_factors(std::string const & factors_str);

}


namespace cutrialdive {

    template <typename PrimeT, typename ExponentT>
    inline
    void sort_factors(std::span<factor<PrimeT, ExponentT>> span)
    {
        std::sort(std::begin(span), std::end(span), [](auto const & x, auto const & y) {
            return x.prime < y.prime;
        });
    }

    template <typename PrimeT, typename ExponentT>
    inline
    std::ostream & operator<<(std::ostream & out, factor<PrimeT, ExponentT> const & factor)
    {
        out << factor.prime;
        if(factor.exponent > ExponentT{1}) {
            out << "^" << factor.exponent;
        }        
        return out;
    }

    template <typename PrimeT, typename ExponentT>
    inline
    std::vector<factor<PrimeT, ExponentT>> parse_factors(std::string const & factors_str)
    {
        using std::operator""sv;
        std::string const factor_regex_str{R"-([1-9][0-9]*(\^[1-9][0-9]*)?)-"};
        std::regex const factors_regex{factor_regex_str + "(,[ ]*" + factor_regex_str + ")*"};
        if(!std::regex_match(factors_str, factors_regex)) {
            std::ostringstream ostr;
            ostr << "Invalid factors list: `" << factors_str << "'";
            throw std::runtime_error{ostr.str()};
        }
        std::vector<factor<PrimeT, ExponentT>> factors;
        constexpr auto delim{","sv};
        for (auto word : std::string_view{factors_str} | std::views::split(',')) {
            auto factor_str = std::string_view(word);
            auto start = factor_str.find(' ');
            if(start != std::string_view::npos) {
                factor_str.remove_prefix(start);
            }
            auto caret_pos = factor_str.find('^');
            PrimeT prime;
            std::istringstream istr_fact{factor_str.substr(0, caret_pos).data()};
            istr_fact >> prime;
            ExponentT exponent{1};
            if(caret_pos != std::string_view::npos) {
                std::istringstream istr_exp{factor_str.substr(caret_pos + 1).data()};
                istr_exp >> exponent;
            }
            factors.emplace_back(prime, exponent);
        }
        return factors;
    }

}
