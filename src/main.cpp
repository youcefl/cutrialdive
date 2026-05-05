/*
* Youcef Lemsafer, 2026.04.21
*/
#include "command_line.hpp"

#include <iostream>
#include <stdexcept>

#include "prp.hpp"
#include "number_sequence.hpp"
#include "builtin_number_sequences.hpp"
#include "trial_factoring.hpp"


namespace {

    using namespace cutrialdive;

    template <typename IndexT>
    int print_value(mode_flag modeFlag, IndexT n, std::ostream& out) {
        return dispatch_mode<IndexT>(modeFlag, [&]<typename Seq>() {
            Seq::print_value(n, out);
            out << std::endl;
            return 0;
        });
    }

    template <typename IndexT>
    int print_expression(mode_flag modeFlag, IndexT n, std::ostream& out) {
        return dispatch_mode<IndexT>(modeFlag, [&]<typename Seq>() {
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
        auto options = *parser.get_options();
        if(options.wants_value) {
            return print_value(options.mode, options.n, std::cout);
        }
        if(options.wants_expression) {
            return print_expression(options.mode, options.n, std::cout);
        }
        if(options.wants_single_prp) {
            run_prp_test(options.mode, options.n, options.factors, options.wants_boosted_factors, std::cout);
            return 0;
        }
//        std::cout << "Trial factoring options: " << *options.tf_options << std::endl;
        trial_factor(options.mode, *options.tf_options);
    } catch(std::exception const & e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    return 0;
}

