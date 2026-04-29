/*
* Youcef Lemsafer, 2026.04.21
*/
#include "siever.hpp"
//#include "hgint.hpp"
#include "smarandache.hpp"

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <cstring>
#include <cstdint>
#include <array>
#include <algorithm>
#include <regex>
#include <charconv>
#include <ranges>


namespace smd = smarandachization;

void output_usage(std::ostream & out)
{
    out << R"-(Usage:
  smarandachization (-p n
                    | --print-expression n 
                    | -s n0 -e n1 --tf-start f0 --tf-end f1
                            [--output <output_path>]
                    | --single-prp n --factors "f_1, ..., f_m"
                    )

  Sm(n) = 123....n (concatenation of the first n integers)

    Print Sm(n)
    -----------
    -p n
            print Sm(n)

    --print-expression n
            prints an expression whose evaluation yields Sm(n)


    Trial factoring
    ---------------
    Use options -s and -e to specify the range of indices [n0, n1[:

    -s n0
            Processing will start at Sm(n0)
    -e n1
            Processing will end at Sm(k) where is the largest integer
            such that n0 <= k < n1
  
    Use options --tf-start --tf-end to specify the range containing
    the prime numbers to consider:

    --tf-start f0
            Lower bound of the set of primes to consider

    --tf-end f1
            Upper bound of the set of primes to consider (f1 is excluded)

    --output <path>
            Causes the output to be written to the given path

    PRP test
    --------
    --single-prp n
            Perform a PRP test on the cofactor of Sm(n).

    If --factors is provided, the factors are divided out before the test.
    Otherwise, the test is performed on Sm(n) itself.

    --factors "f1, f2, ..., fm"
            A comma-separated list of known factors of Sm(n).

            Each factor can be expressed as:
            - A simple prime: "p"
            - A prime with exponent: "p^e"

            Example: "2^5, 3^2, 223, 1279"
            Spaces after commas are optional.

    --boost-factors
            When enabled, for each factor p, the program automatically
            determines the largest exponent e such that p^e divides Sm(n).
            This ensures that the cofactor used for the PRP test is not
            divisible by any of the given primes.

)-";
}

template <typename T>
T expect_value(char const* optName, char** parg)
{
    if(!*parg) {
        std::string msg{"Missing value for option "};
        msg = msg + optName;
        throw std::runtime_error(msg.c_str());
    }
    if constexpr (std::is_same_v<std::string, T>) {
        return std::string{*parg};
    }
    T val;
    std::istringstream{*parg} >> val;
    
    return val;
}

template <typename T, typename Action>
bool eat_valued_option(char const* flag, char**& argv, Action act)
{
    if(!*argv) {
        throw std::runtime_error{std::string{"Invalid command line, missing option: "} + flag};
    }
    if(!strcmp(*argv, flag)) {
        auto val = expect_value<T>(flag, ++argv);
        act(std::move(val));
        return true;
    }
    return false;
}

std::vector<smd::factor<uint64_t>> parse_factors(std::string const & factors_str)
{
    using std::operator""sv;
    std::string const factor_regex_str{R"-([1-9][0-9]*(\^[1-9][0-9]*)?)-"};
    std::regex const factors_regex{factor_regex_str + "(,[ ]*" + factor_regex_str + ")*"};
    if(!std::regex_match(factors_str, factors_regex)) {
        throw std::runtime_error{"Invalid factors list"};
    }
    std::vector<smd::factor<uint64_t>> factors;
    constexpr auto delim{","sv};
    for (auto word : std::string_view{factors_str} | std::views::split(',')) {
        auto factor_str = std::string_view(word);
        auto start = factor_str.find(' ');
        if(start != std::string_view::npos) {
            factor_str.remove_prefix(start);
        }
        auto caret_pos = factor_str.find('^');
        uint64_t prime;
        std::istringstream istr_fact{factor_str.substr(0, caret_pos).data()};
        istr_fact >> prime;
        uint64_t exponent{1};
        if(caret_pos != std::string_view::npos) {
            std::istringstream istr_exp{factor_str.substr(caret_pos + 1).data()};
            istr_exp >> exponent;
        }
        factors.emplace_back(prime, exponent);
    }
    return factors;
}


int main(int argc, char** argv)
{
    try {
        if(argc <= 2) {
            output_usage(std::cout);
            return 1;
        }
        bool haveStart{}, haveEnd{}, haveTfStart{}, haveTfEnd{}, haveOutputPath;
        uint64_t startIndex{}, endIndex{},
                tfStart{}, tfEnd{};
        std::string output_path;
        for(auto ppargv = argv + 1; *ppargv; ++ppargv) {
            if(eat_valued_option<uint64_t>("-p", ppargv, [](auto index) {
                smd::output_smarandache(std::cout, index) << std::endl;
            })) {
                return 0;
            }
            if(eat_valued_option<uint64_t>("--print-expression", ppargv, [](auto index){
                smd::output_smarandache_expression(std::cout, index) << std::endl;
            })){
                return 0;
            }
            uint64_t index{};
            if(eat_valued_option<uint64_t>("--single-prp", ppargv, [&index](auto idx) { index = idx; })) {
                std::vector<smd::factor<uint64_t>> factors;
                ++ppargv; // go past value of --single-prp
                if(eat_valued_option<std::string>("--factors", ppargv, [&factors](auto factors_str) {
                    factors = parse_factors(factors_str);
                })) {
                    ++ppargv;
                }
                auto haveToBoostFactors = false;
                if(*ppargv && !strcmp(*ppargv, "--boost-factors")) {
                    haveToBoostFactors = true;
                }
                smd::run_prp_test(index, factors, haveToBoostFactors);
                if(haveToBoostFactors) {
                    std::cout << "Factors used: ";
                    char const * sep = "";
                    std::for_each(std::begin(factors), std::end(factors), [&sep](auto factor){
                        std::cout << sep << factor.value;
                        if(factor.exponent > 1) {
                            std::cout << "^" << factor.exponent;
                        }
                        sep = ", ";
                    });
                    std::cout << (*sep ? "" : "<none>") << std::endl;
                }
                return 0;
            }
            haveStart = haveStart || eat_valued_option<uint64_t>("-s", ppargv, 
                                        [&startIndex](auto start){ startIndex = start; });
            haveEnd =  haveEnd || eat_valued_option<uint64_t>("-e", ppargv,
                                        [&endIndex](auto end){ endIndex = end; });
            haveTfStart = haveTfStart || eat_valued_option<uint64_t>("--tf-start", ppargv,
                                        [&tfStart](auto start){ tfStart = start; });
            haveTfEnd = haveTfEnd || eat_valued_option<uint64_t>("--tf-end", ppargv,
                                        [&tfEnd](auto end){ tfEnd = end; });
            haveOutputPath = eat_valued_option<std::string>("--output", ppargv,
                                        [&output_path](auto path){ output_path = path; });
        }
        auto flags = std::array{haveStart, haveEnd, haveTfStart, haveTfEnd};
        if(std::any_of(std::begin(flags), std::end(flags), [](auto val){ return !val; })
            || (haveOutputPath && output_path.empty())) {
            std::cerr << "Invalid command line" << std::endl;
            return 1;
        }
        std::vector<smd::result<uint64_t>> results{};
        smd::time("Overall time: ", [&](){
            smd::trial_factor(smd::trial_factoring_options{startIndex, endIndex,
                tfStart, tfEnd,
                haveOutputPath ? std::filesystem::path{output_path} : std::optional<std::filesystem::path>{}
            });
        });
    } catch(std::exception const & e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

