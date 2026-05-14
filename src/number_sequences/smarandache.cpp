/*
 * Creation date: 2026.04.21
 * Created by Youcef Lemsafer
 */
#include "smarandache.hpp"

#include <cstdint>
#include <algorithm>
#include <ostream>
#include <array>
#include <functional>

#include "hgint.hpp"


namespace cutrialdive
{
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

    /// @brief Returns Sm(n) = 1234567...n
    /// @tparam NumberT number type
    /// @param n 
    template <typename T>
    T smarandache(uint64_t n)
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

    /// @brief Outputs Sm(n)
    /// @param[in] n indice
    /// @param[in,out] out the output stream
    /// @return output stream
    std::ostream &output_smarandache(uint64_t n, std::ostream &out)
    {
        if(!n) {
            out << "0";
            return out;
        }
        for(decltype(n) i = 1; i <= n; ++i) {
            out << i;
        }
        return out;
    }

    /// @brief Outputs Sm(n) as an expression
    /// @param[in/out] out the output stream
    /// @param[in] n indice
    /// @return output stream
    std::ostream &output_smarandache_expression(uint64_t n, std::ostream &out)
    {
        if(n <= 1) {
            return out << n;
        }
        auto log10n = details::floor_log10(n);
        static const auto c3 = std::string{"1490840987654358*10^180-10990220"};
        static const auto c4 = std::string{"10^2701*("} + c3 + ")*10201-1109890222000";
        static const auto c5 = std::string{"10^36000*("} + c4 + ")*123454321-123443211358*33300^2";
        static const std::array<std::function<std::string(std::string)>, 5> sm = {
            [](auto x) { return std::string{"(10^("} + x + "+1)"
                                        + "-(9*" + x + "+10))/81"; }
          , [](auto x) { return std::string{"(10^(2*"} + x + "-18)*(123456789*99^2+991)"
                            "-(99*" + x + "+100))/9801"; }
          , [](auto x) { return std::string{"(("} + c3 + ")*10^(3*" + x + "-296)"
                            "-(120879*" + x + "+121000))/120758121"; }
          , [](auto x) { return std::string{"(("} + c4 + ")*10^(4*(" + x + "-999))"
                            "-12321*(9999*" + x + "+10000))/1231853592321"; }
          
          , [](auto x) { return std::string{"(("} + c5 + ")*10^(5*(" + x + "-9999))"
                                              "-12321*1234321*(99999*" + x + "+100000))"
                                            "/(12321*1234321*99999^2)"; }
        };
        if(log10n < sm.size()) {
            return out << sm[log10n](std::to_string(n));
        }
        throw std::runtime_error{"smarandache(n) for n >= 10^5 not supported yet"};
    }

    char const*
    number_sequence<num_seq_id::smarandache>::short_name()
    {
        return "Sm";
    }

    typename number_sequence<num_seq_id::smarandache>::value_type
    number_sequence<num_seq_id::smarandache>::value(index_type n)
    {
        return smarandache<value_type>(n);
    }

    void
    number_sequence<num_seq_id::smarandache>::print_value(index_type n, std::ostream & out)
    {
        output_smarandache(n, out);
    }

    void
    number_sequence<num_seq_id::smarandache>::print_expression(index_type n, std::ostream & out)
    {
        output_smarandache_expression(n, out);
    }
}
