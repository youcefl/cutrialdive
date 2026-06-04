/*
* Created on 2026.05.04
* Copyright (c) Youcef Lemsafer
*/
#include "command_line.hpp"

#include <string.h>
#include <string>
#include <sstream>

#include "common_defines.h"
#include "factor.hpp"

namespace {
    using namespace cutrialdive;

    template <typename T, typename ExceptionBuilder>
    inline
    T expect_value(ExceptionBuilder buildException, char const* optName, char**& parg)
    {
        if(!*parg) {
            throw buildException(std::string{"Missing value for option "} + optName);
        }
        if constexpr (std::is_same_v<std::string, T>) {
            return std::string{*parg++};
        }
        T val;
        std::istringstream{*parg} >> val;
        ++parg;
        return val;
    }

    template <typename T, typename Action, typename ExceptionBuilder>
    inline
    bool eat_valued_option(char const* flag, char**& argv, ExceptionBuilder buildException, Action act)
    {
        if(!*argv) {
            return false;
        }
        if(!strcmp(*argv, flag)) {
            ++argv;
            auto val = expect_value<T>(buildException, flag, argv);
            act(std::move(val));
            return true;
        }
        return false;
    }

    template <typename ExceptionBuilder>
    inline
    void parse_num_seq(command_line_options & options, char**& argv, ExceptionBuilder buildException)
    {
        if(!eat_valued_option<std::string>("--num-seq", argv, buildException, [&](auto numSeqSpec) {
                try {
                    auto commaPos = numSeqSpec.find(',');
                    if(commaPos != std::string::npos) {
                        auto numSeqIdAsStr = numSeqSpec.substr(0, commaPos);
                        options.num_seq_id_value = num_seq_id_from_string(numSeqIdAsStr);
                        options.num_seq_params = numSeqSpec.substr(commaPos + 1);
                        if(options.num_seq_params.empty()) {
                            throw buildException("Invalid value for option --num-seq");
                        }
                    } else {
                        options.num_seq_id_value = num_seq_id_from_string(numSeqSpec);
                    }
                } catch (std::exception const & ex) {
                    throw buildException(ex.what());
                }
            })) {
           throw buildException("Invalid command line, option --num-seq must come first");
        }
    }
}

namespace cutrialdive {

    command_line_parse_exception::command_line_parse_exception(std::string const & msg)
        : std::runtime_error(msg)
    {}

    command_line_parse_exception::command_line_parse_exception(char const* msg)
        : std::runtime_error(msg)
    {}


    command_line_parser::command_line_parser(int argc, char** argv)
        : wants_help_{}
        , wants_usage_{}
        , wants_version_{}
        , wants_autotests_{}
    {
        if(argc == 2) {
            if(!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
                wants_help_ = true;
                return;
            }
            if(!strcmp(argv[1], "--version")) {
                wants_version_ = true;
                return;
            }
            if(!strcmp(argv[1], "--autotest")) {
                wants_autotests_ = true;
                return;
            }
        }
        if(argc <= 4) {
            wants_usage_ = true;
            return;
        }
        auto ppargv = argv + 1;

        trial_factoring_options tfOptions{};
        options_ = command_line_options{};
        auto buildException = [](auto msg){ return command_line_parse_exception{msg}; };
        parse_num_seq(*options_, ppargv, buildException);

        bool haveStart{}, haveEnd{}, haveTfBits{}, haveTfStart{}, haveTfEnd{};
        for(; *ppargv; ++ppargv) {
            if(eat_valued_option<uint64_t>("--print", ppargv, buildException, [this](auto n) {
                options_->wants_value = true;
                options_->n = n;
            })) {
                return;
            }
            if(eat_valued_option<uint64_t>("--print-expression", ppargv, buildException, [this](auto n){
                options_->wants_expression = true;
                options_->n = n;
            })){
                return;
            }
#ifdef CUTRIALDIVE_ENABLE_PRP
            if(eat_valued_option<uint64_t>("--single-prp", ppargv, buildException, [this](auto n) {
                options_->wants_single_prp = true;
                options_->n = n;
            })) {
                std::vector<factor<uint64_t, uint32_t>> factors;
                eat_valued_option<std::string>("--factors", ppargv, buildException, [&](auto factors_str) {
                    options_->factors = parse_factors<uint64_t, uint32_t>(factors_str);
                });
                if(*ppargv && !strcmp(*ppargv, "--no-boost-factors")) {
                    options_->wants_boosted_factors = false;
                }
                return;
            }
#endif // CUTRIALDIVE_ENABLE_PRP
            haveStart = haveStart || eat_valued_option<uint64_t>("-s", ppargv, buildException,
                                        [&tfOptions](auto n0){ tfOptions.n0 = n0; });
            haveEnd =  haveEnd || eat_valued_option<uint64_t>("-e", ppargv, buildException,
                                        [&tfOptions](auto n1){ tfOptions.n1 = n1; });
            haveTfBits = haveTfBits || eat_valued_option<uint16_t>("--tf-bits", ppargv, buildException,
                [&](auto tfBits){
                    tfOptions.f0 = 2;
                    if(tfBits > 64) {
                        throw command_line_parse_exception{
                            "The value of option --tf-bits must be an integer between 0 and 64"};
                    }
                    tfOptions.f1 = (tfBits < 64) ? (uint64_t{1} << tfBits) : ~uint64_t{};
                });
            haveTfStart = haveTfStart || eat_valued_option<uint64_t>("--tf-start", ppargv, buildException,
                                        [&tfOptions](auto start){ tfOptions.f0 = start; });
            haveTfEnd = haveTfEnd || eat_valued_option<uint64_t>("--tf-end", ppargv, buildException,
                                        [&tfOptions](auto end){ tfOptions.f1 = end; });
            if(haveTfBits && (haveTfStart || haveTfEnd)) {
                throw command_line_parse_exception{
                            "Mixing --tf-bits with --tf-start and --tf-end is not allowed"};
            }
            eat_valued_option<std::string>("--output", ppargv, buildException,
                        [&tfOptions](auto path){ tfOptions.output_path = path; });
        }
        auto mandatoryFlags = std::array{haveStart, haveTfStart || haveTfBits, haveTfEnd || haveTfBits};
        if(std::any_of(std::begin(mandatoryFlags), std::end(mandatoryFlags), [](auto val){ return !val; })) {
            throw command_line_parse_exception{"Missing mandatory option"};
        }
        if(!haveEnd) {
            tfOptions.n1 = tfOptions.n0 + 1;
        }
        if (tfOptions.output_path && tfOptions.output_path->empty()) {
            throw command_line_parse_exception{"Invalid output path"};
        }
        options_->tf_options = std::move(tfOptions);
    }

    void print_version(std::ostream & out)
    {
        out << "cutrialdive version " << CUTRIALDIVE_VERSION << std::endl;
    }

    void print_usage_only(std::ostream & out)
    {
        out << R"-(
  cutrialdive ( --version
              | ( -h | --help )
              | ( --num-seq <NUM_SEQ_ID>
                        ( --print n
                        | --print-expression n 
                        | -s n0 [-e n1] ( --tf-bits bits_count 
                                        | --tf-start f0 --tf-end f1 )
                          [--output <output_path>]
)-"
#ifdef CUTRIALDIVE_ENABLE_PRP
R"-(                        | --single-prp n [--factors "f_1, ..., f_m"
                                                [--no-boost-factors]]
)-"
#endif // CUTRIALDIVE_ENABLE_PRP
R"-(                        )
                )
              )
)-"
        ;
    }

    void print_usage(std::ostream & out)
    {
        print_usage_only(out);
        out << R"-(

  Type 'cutrialdive --help' for more.
)-";
    }

    void print_help(std::ostream & out)
    {
        print_usage_only(out);
        out << R"-(

    --help
            prints this help and exits.
    --version
            prints the version number and exits.
    
    Number sequence selection
    -------------------------
    --num-seq <NUM_SEQ_ID>
            Selects a number sequence to work on, possible values of
            <NUM_SEQ_ID> are: smarandache[,base], mersenne.
            smarandache,2 means Smarandache base 2, default is base 10.

    Print number
    ------------
    --print n
            prints S(n) in base 10

    --print-expression n
            prints an expression whose evaluation yields S(n). For example
            "--num-seq mersenne --print-expression 7" will output 2^7-1.

    Trial factoring
    ---------------
    Use options -s and -e to specify the range of indices [n0, n1[:

    -s n0
            Processing will start with S(n0)
    -e n1
            Processing will end with S(k) where k is the largest integer such
            that n0 <= k < n1.
            Default value is n0 + 1
  
    Use option --tf-bits or --tf-start and --tf-end to specify the range
    to sieve for prime numbers to consider:

    --tf-bits  bits_count
            Use the primes in [0, 2^bits_count[

    --tf-start f0
            Lower bound of the set of primes to consider

    --tf-end f1
            Upper bound of the set of primes to consider (f1 is excluded)

    --output <path>
            Causes the output to be written to the given path
)-"
#ifdef CUTRIALDIVE_ENABLE_PRP
R"-(
    PRP test
    --------
    --single-prp n
            Perform a PRP test on the cofactor of S(n).

    If --factors is provided, the factors are divided out before the test.
    Otherwise, the test is performed on S(n) itself.

    --factors "f1, f2, ..., fm"
            A comma-separated list of known factors of F(n).

            Each factor can be expressed as:
            - A simple prime: "p"
            - A prime with exponent: "p^e"

            Example: "2^5, 3^2, 223, 1279"
            Spaces after commas are optional.

    --no-boost-factors
            Consider the factors has already having the maximal exponents.
            By default, for each factor p, the program automatically
            determines the largest exponent e such that p^e divides S(n).
            This ensures that the cofactor used for the PRP test is not
            divisible by any of the given primes.
            This option disables that behavior.
)-"
#endif // CUTRIALDIVE_ENABLE_PRP
R"-(
)-"        ;
    }

}

