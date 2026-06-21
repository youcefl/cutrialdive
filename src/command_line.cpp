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
    T expect_value(ExceptionBuilder buildException, char const* optName, char** parg)
    {
        if(!*parg) {
            throw buildException(std::string{"Missing value for option "} + optName);
        }
        if constexpr (std::is_same_v<std::string, T>) {
            return std::string{*parg};
        }
        T val;
        std::istringstream{*parg} >> val;
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

    inline
    bool parse_no_progress(command_line_options & options, char**& argv)
    {
        if(!strcmp(*argv, "--no-progress")) {
            options.is_progress_enabled = false;
            return true;
        }
        return false;
    }
    template <typename ValueT, typename DestT, typename ExceptionBuilder>
    inline
    bool parse_valued_option(char const * optName, DestT & dest, char**& argv, ExceptionBuilder buildException)
    {
        return eat_valued_option<ValueT>(optName, argv, buildException, [&](auto value){
            dest = value;
        });
    }

    template <typename ExceptionBuilder>
    inline
    bool parse_verbosity(command_line_options & options, char**& argv, ExceptionBuilder buildException)
    {
        return eat_valued_option<uint16_t>("--verbosity", argv, buildException, [&](auto verbosityLevel){
                if(!is_valid_verbosity(verbosityLevel)) {
                    throw buildException("Invalid value for option --verbosity");
                }
                options.verbosity_level = static_cast<verbosity>(verbosityLevel);
            });
    }

    template <typename ExceptionBuilder>
    inline
    bool parse_threads(command_line_options & options, char**& argv, ExceptionBuilder buildException)
    {
        return eat_valued_option<uint32_t>("--threads", argv, buildException, [&](auto threadsCount){
                options.threads_count = threadsCount;
            });
    }

    template <typename ExceptionBuilder>
    inline
    bool parse_segment_length(command_line_options & options, char**& argv, ExceptionBuilder buildException)
    {
        if(eat_valued_option<uint64_t>("--segment-length", argv, buildException, [&](auto segmentLength){
            if(options.segment_length.has_value()) {
                throw buildException(
                        "Invalid command line, --segment-length is not allowed when --segment-bits is present");
            }
            options.segment_length = segmentLength;
        })) {
            return true;
        }
        if(eat_valued_option<int32_t>("--segment-bits", argv, buildException, [&](auto segmentBits){
            if(options.segment_length.has_value()) {
                throw buildException(
                        "Invalid command line, --segment-bits is not allowed when --segment-length is present");
            }
            if(segmentBits < 0 || segmentBits >= 64) {
                throw buildException(
                        "value of option --segment-bits must be an integer in range [0, 63]");
            }
            options.segment_length = uint64_t{1} << segmentBits;
        })) {
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
                    options.seq_id = num_seq_id_from_string(numSeqIdAsStr);
                    options.seq_params = numSeqSpec.substr(commaPos + 1);
                    if(options.seq_params.empty()) {
                        throw buildException("Invalid value for option --num-seq");
                    }
                } else {
                    options.seq_id = num_seq_id_from_string(numSeqSpec);
                }
            } catch (std::exception const & ex) {
                throw buildException(ex.what());
            }
            ++argv;
        })) {
           throw buildException("Invalid command line, option --num-seq must come first");
        }
    }

#ifdef CUTRIALDIVE_ENABLE_GPU
    template <typename ExceptionBuilder>
    inline
    bool parse_device_options(command_line_options & options, char**& argv, ExceptionBuilder buildException)
    {
        if(parse_valued_option<int32_t>("-g", options.device_id, argv, buildException)) {
            return true;
        }
        if(parse_valued_option<int32_t>("--grid-size", options.grid_size, argv, buildException)) {
            return true;
        }
        if(parse_valued_option<int32_t>("--block-size", options.block_size, argv, buildException)) {
            return true;
        }
        return false;
    }
#endif

    template <typename ExceptionBuilder>
    inline
    void parse_resume(command_line_options & options, char**& argv, ExceptionBuilder buildException)
    {
        options.is_resuming = false;
        bool haveChkpntPeriod = false;
        for(; *argv; ++argv) {
            if(eat_valued_option<std::string>("--resume", argv, buildException, [&](auto chkpntPath){
                options.is_resuming = true;
                options.checkpoint_path = chkpntPath;
            })) {
                continue;
            }
#ifdef CUTRIALDIVE_ENABLE_GPU
            if(parse_device_options(options, argv, buildException)) {
                continue;
            }
#endif
            if(eat_valued_option<uint32_t>("--checkpoint-period", argv, buildException, [&](auto chkpntPeriod){
                options.checkpoint_period = chkpntPeriod;
                haveChkpntPeriod = true;
            })) {
                continue;
            }
            if(parse_verbosity(options, argv, buildException)) {
                continue;
            }
            if(parse_valued_option<uint32_t>("--threads", options.threads_count, argv, buildException)) {
                continue;
            }
            if(parse_no_progress(options, argv)) {
                continue;
            }
            if(parse_segment_length(options, argv, buildException)) {
                continue;
            }
            if(!options.is_resuming && !haveChkpntPeriod) {
                return;
            }
        }
        if(haveChkpntPeriod && !options.is_resuming) {
            throw buildException("Invalid command line, option --checkpoint-period is unexpected");
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
        if(argc == 1) {
            wants_usage_ = true;
            return;
        }
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

        auto ppargv = argv + 1;
        auto buildException = [](auto msg){ return command_line_parse_exception{msg}; };
        options_ = command_line_options{};

        parse_resume(*options_, ppargv, buildException);
        if(options_->is_resuming) {
            return;
        }

        if(parse_verbosity(*options_, ppargv, buildException)) {
            ++ppargv;
        }
        parse_num_seq(*options_, ppargv, buildException);

        trial_factoring_options tfOptions{};
        bool haveStart{}, haveEnd{}, haveTfBits{}, haveTfStart{}, haveTfEnd{};

        for(; *ppargv; ++ppargv) {
            if(parse_verbosity(*options_, ppargv, buildException)) {
                continue;
            }
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
#ifdef CUTRIALDIVE_ENABLE_GPU
            if(parse_device_options(*options_, ppargv, buildException)) {
                continue;
            }
#endif
            if(parse_valued_option<uint32_t>("--threads", options_->threads_count, ppargv, buildException)) {
                continue;
            }
            if(parse_valued_option<uint64_t>("-s", tfOptions.n0, ppargv, buildException)) {
                haveStart = true;
                continue;
            }
            if(parse_valued_option<uint64_t>("-e", tfOptions.n1, ppargv, buildException)) {
                haveEnd = true;
                continue;
            }
            if(eat_valued_option<uint16_t>("--tf-bits", ppargv, buildException,
                [&](auto tfBits){
                    tfOptions.f0 = 2;
                    if(tfBits > 64) {
                        throw command_line_parse_exception{
                            "The value of option --tf-bits must be an integer between 0 and 64"};
                    }
                    tfOptions.f1 = (tfBits < 64) ? (uint64_t{1} << tfBits) : ~uint64_t{};
                })){
                    haveTfBits = true;
                    continue;
            }
            if(parse_valued_option<uint64_t>("--tf-start", tfOptions.f0, ppargv, buildException)) {
                haveTfStart = true;
                continue;
            }
            if(parse_valued_option<uint64_t>("--tf-end", tfOptions.f1, ppargv, buildException)) {
                haveTfEnd = true;
                continue;
            }
            if(haveTfBits && (haveTfStart || haveTfEnd)) {
                throw command_line_parse_exception{
                            "Mixing --tf-bits with --tf-start and --tf-end is not allowed"};
            }
            if(parse_no_progress(*options_, ppargv)) {
                continue;
            }
            if(parse_valued_option<std::string>("--output", tfOptions.output_path, ppargv, buildException)) {
                continue;
            }
            if(parse_valued_option<uint32_t>("--checkpoint-period", options_->checkpoint_period, ppargv, buildException)) {
                continue;
            }
            if(parse_segment_length(*options_, ppargv, buildException)) {
                continue;
            }
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
                        | <trial factoring options>
)-"
#ifdef CUTRIALDIVE_ENABLE_PRP
R"-(                        | <prp test options>
)-"
#endif // CUTRIALDIVE_ENABLE_PRP
R"-(                        )
                )
              | --resume <checkpoint_path> <resume options>
              )
              [<global options>]
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
    
    Global options
    --------------
    --threads <N>
            Use up to N concurrent threads.

    --verbosity <V>
            Set log verbosity level to V where V is an integer in range [0, 4]
            0: silent (no messages)
            1: normal
            2: verbose
            3: debug
            4: diagnostic
    
    Number sequence selection
    -------------------------
    --num-seq <NUM_SEQ_ID>
            Selects a number sequence to work on, possible values of
            <NUM_SEQ_ID> are: smarandache[,base], mersenne.
            smarandache,2 means Smarandache base 2, default is base 10.

    Print number
    ------------
    --print n
            prints S(n) in base 10.

    --print-expression n
            prints an expression whose evaluation yields S(n). For example
            "--num-seq mersenne --print-expression 7" will output 2^7-1.

    Trial factoring
    ---------------
    If no output path is specified (i.e. option --output is missing) then
    the results are written to the console and checkpointing is deactivated.
    If an output path P is specified then P.chkpnt is reserved for the
    checkpoint file. The checkpoint file is deleted upon succesful completion
    of the trial factoring.

    Use options -s and -e to specify the range of indices [n0, n1[:

    -s n0
            Processing will start with S(n0).
    -e n1
            Processing will end with S(k) where k is the largest integer such
            that n0 <= k < n1.
            Default value is n0 + 1
  
    Use option --tf-bits or --tf-start and --tf-end to specify the range
    to sieve for prime numbers to consider:

    --tf-bits bits_count
            Use the primes in [0, 2^bits_count[.

    --tf-start f0
            Lower bound of the set of primes to consider.

    --tf-end f1
            Upper bound of the set of primes to consider (f1 is excluded).

    --no-progress
            Suppresses progress display.

    --output <path>
            Causes the output to be written to the given path.

    --checkpoint-period T
            Minimum duration in seconds between two checkpoint file savings.

    --segment-bits <L>
            Use 2^L as the length of the sieve segment, not compatible
            with --segment-length
    --segment-length <N>
            Use N as the length of the sieve segment, not compatible
            with --segment-bits.
            The default sieve segment length is 2^26.
)-"
#ifdef CUTRIALDIVE_ENABLE_GPU
R"-(
    -g <device_id>
            Id of the device to run the trial factoring on.

    --grid-size N
            Override default grid size (blocks per kernel launch) with N.

    --block-size N
            Override default block size (threads per block) with N.
)-"
#endif
R"-(
    Resuming a trial factoring
    --------------------------
    An interrupted trial factoring can be resumed with option --resume.
    In the resume mode some of the options above can be specified/altered.

    --resume <checkpoint file path>
            Resume an interrupted trial factoring.
    
    The options available in resume mode are the following: --threads,
    --no-progress, --checkpoint-period)-"
#ifdef CUTRIALDIVE_ENABLE_GPU
                                  R"-(, -g, --grid-size, --block-size
)-"
#endif
R"-(        and any of the global options.
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
)-";
    }

}

