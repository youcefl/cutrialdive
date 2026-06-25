/*
* MIT License
* Created on 2026.04.21
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
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
#include "logger.hpp"

namespace {
    using namespace cutrialdive;

    template <typename IndexT>
    int print_value(num_seq_spec numSeqSpec, IndexT n, std::ostream& out)
    {
        return dispatch_num_seq<IndexT>(numSeqSpec, [&]<typename Seq, typename... Args>(Args&&... args) {
            if constexpr(HasValuePrinter<Seq>) {
                get_math_sequence(Seq{std::forward<Args>(args)...}).print_value(n, out) << std::endl;
            } else {
                out << "Value not available" << std::endl;
            }
            return 0;
        });
    }

    template <typename IndexT>
    int print_expression(num_seq_spec numSeqSpec, IndexT n, std::ostream& out)
    {
        return dispatch_num_seq<IndexT>(numSeqSpec, [&]<typename Seq, typename... Args>(Args&&... args) {
            if constexpr(HasExpressionPrinter<Seq>) {
                get_math_sequence(Seq{std::forward<Args>(args)...}).print_expression(n, out) << std::endl;
            } else {
                out << "Expression not available" << std::endl;
            }
            return 0;
        });
    }

    runtime_options get_runtime_options(command_line_options const & cmdOpts)
    {
        auto rtOpts = runtime_options::default_options();
        if(cmdOpts.threads_count.has_value()) {
            rtOpts.threads_count = *cmdOpts.threads_count;
        }
        if(cmdOpts.verbosity_level.has_value()) {
            rtOpts.verbosity_level = *cmdOpts.verbosity_level;
        }
        return rtOpts;
    }

    tf_runtime_options get_tf_runtime_options(command_line_options const & cmdOpts)
    {
        auto tfRuntimeOpts = tf_runtime_options::default_options();
        tfRuntimeOpts.common_options = get_runtime_options(cmdOpts);
        if(cmdOpts.checkpoint_period) {
            tfRuntimeOpts.checkpoint_period = std::chrono::seconds{*cmdOpts.checkpoint_period};
        }
        if(cmdOpts.is_progress_enabled.has_value()) {
            tfRuntimeOpts.is_progress_enabled = *cmdOpts.is_progress_enabled;
        }
        if(cmdOpts.segment_length.has_value()) {
            tfRuntimeOpts.segment_length = *cmdOpts.segment_length;
        }
#ifdef CUTRIALDIVE_ENABLE_GPU
        if(cmdOpts.device_id.has_value()) {
            tfRuntimeOpts.device_id = *cmdOpts.device_id;
        }
        if(cmdOpts.grid_size.has_value()) {
            tfRuntimeOpts.grid_size = *cmdOpts.grid_size;
        }
        if(cmdOpts.block_size.has_value()) {
            tfRuntimeOpts.block_size = *cmdOpts.block_size;
        }
#endif
        return tfRuntimeOpts;
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
        auto cmdOptions = *parser.get_options();
        auto rtOpts = get_runtime_options(cmdOptions);
        log().set_verbosity(rtOpts.verbosity_level);
        log().set_stream(std::cout);
        if(cmdOptions.is_resuming) {
            ctd::resume_trial_factoring(cmdOptions.checkpoint_path.value(), get_tf_runtime_options(cmdOptions), std::cout);
            return 0;
        }
        auto numSeqSpec = num_seq_spec{cmdOptions.seq_id, cmdOptions.seq_params};
        if(cmdOptions.wants_value) {
            return print_value(numSeqSpec, cmdOptions.n, std::cout);
        }
        if(cmdOptions.wants_expression) {
            return print_expression(numSeqSpec, cmdOptions.n, std::cout);
        }
#ifdef CUTRIALDIVE_ENABLE_PRP
        if(cmdOptions.wants_prp) {
            ctd::run_prp_test(
                numSeqSpec,
                cmdOptions.n,
                cmdOptions.factors,
                cmdOptions.wants_boosted_factors,
                std::cout
            );
            return 0;
        }
#endif
        ctd::trial_factor(numSeqSpec, *cmdOptions.tf_options, get_tf_runtime_options(cmdOptions), std::cout);
    } catch(std::exception const & e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    return 0;
}

