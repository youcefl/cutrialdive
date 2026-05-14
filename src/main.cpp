/*
* Youcef Lemsafer, 2026.04.21
*/
#include "command_line.hpp"

#include <iostream>
#include <stdexcept>

#ifdef CUTRIALDIVE_ENABLE_PRP
#include "prp.hpp"
#endif
#include "builtin_number_sequences.hpp"
#include "num_seq_dispatch.hpp"
#include "trial_factoring.hpp"
#include "autotests.hpp"

namespace {

    using namespace cutrialdive;

    template <typename IndexT>
    int print_value(num_seq_id numSeqId, IndexT n, std::ostream& out) {
        return dispatch_num_seq<IndexT>(numSeqId, [&]<typename Seq>() {
            Seq::print_value(n, out);
            out << std::endl;
            return 0;
        });
    }

    template <typename IndexT>
    int print_expression(num_seq_id numSeqId, IndexT n, std::ostream& out) {
        return dispatch_num_seq<IndexT>(numSeqId, [&]<typename Seq>() {
            Seq::print_expression(n, out);
            out << std::endl;
            return 0;
        });
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
        if(options.wants_value) {
            return print_value(options.num_seq_id_value, options.n, std::cout);
        }
        if(options.wants_expression) {
            return print_expression(options.num_seq_id_value, options.n, std::cout);
        }
#ifdef CUTRIALDIVE_ENABLE_PRP
        if(options.wants_single_prp) {
            run_prp_test(options.num_seq_id_value, options.n, options.factors, options.wants_boosted_factors, std::cout);
            return 0;
        }
#endif
        trial_factor(options.num_seq_id_value, *options.tf_options);
    } catch(std::exception const & e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    return 0;
}

