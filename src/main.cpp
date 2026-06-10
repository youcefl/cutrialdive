/*
* Youcef Lemsafer, 2026.04.21
*/
#include "command_line.hpp"

#include <iostream>
#include <stdexcept>

#ifdef CUTRIALDIVE_ENABLE_PRP
#include "prp.hpp"
#endif
#include "number_sequence.hpp"
#include "builtin_number_sequences.hpp"
#include "num_seq_spec.hpp"
#include "num_seq_dispatch.hpp"
#include "trial_factoring.hpp"
#include "autotests.hpp"

namespace {
    using namespace cutrialdive;

    template <typename IndexT>
    int print_value(num_seq_spec numSeqSpec, IndexT n, std::ostream& out)
    {
        return dispatch_num_seq<IndexT>(numSeqSpec, [&]<typename Seq>() {
            if constexpr(HasValuePrinter<Seq>) {
                Seq{}.print_value(n, out) << std::endl;
            } else {
                out << "Value not available" << std::endl;
            }
            return 0;
        });
    }

    template <typename IndexT>
    int print_expression(num_seq_spec numSeqSpec, IndexT n, std::ostream& out)
    {
        return dispatch_num_seq<IndexT>(numSeqSpec, [&]<typename Seq>() {
            if constexpr(HasExpressionPrinter<Seq>) {
                Seq{}.print_expression(n, out) << std::endl;
            } else {
                out << "Expression not available" << std::endl;
            }
            return 0;
        });
    }

    tf_runtime_options get_runtime_options(command_line_options const & cmdOpts)
    {
        auto runtimeOpts = tf_runtime_options::default_options();
        if(cmdOpts.checkpoint_period) {
            runtimeOpts.checkpoint_period = std::chrono::seconds{*cmdOpts.checkpoint_period};
        }
        return runtimeOpts;
    }
}

namespace ctd = cutrialdive;

int main(int argc, char** argv)
{
    try {
        ctd::command_line_parser parser{argc, argv};
        if(parser.is_version_requested()) {
            print_version(std::cout);
            return 0;
        }
        if(parser.is_usage_requested()) {
            print_usage(std::cerr);
            return 1;
        }
        if(parser.is_help_requested()) {
            print_help(std::cout);
            return 0;
        }
        if(parser.are_autotests_requested()) {
            return ctd::autotest();
        }
        auto options = *parser.get_options();
        if(options.is_resuming) {
            ctd::resume_trial_factoring(options.checkpoint_path.value(), get_runtime_options(options), std::cout);
        }
        auto numSeqSpec = num_seq_spec{options.seq_id, options.seq_params};
        if(options.wants_value) {
            return print_value(numSeqSpec, options.n, std::cout);
        }
        if(options.wants_expression) {
            return print_expression(numSeqSpec, options.n, std::cout);
        }
#ifdef CUTRIALDIVE_ENABLE_PRP
        if(options.wants_single_prp) {
            ctd::run_prp_test(numSeqSpec,
                options.n,
                options.factors,
                options.wants_boosted_factors,
                std::cout
            );
            return 0;
        }
#endif
        ctd::trial_factor(numSeqSpec, *options.tf_options, get_runtime_options(options), std::cout);
    } catch(std::exception const & e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    return 0;
}

