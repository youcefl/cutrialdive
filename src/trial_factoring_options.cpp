/*
* Created on 2026.05.05
* Copyright (c) Youcef Lemsafer
*/
#include "trial_factoring_options.hpp"

namespace cutrialdive {

    std::ostream & operator<<(std::ostream & out, trial_factoring_options const & options)
    {
        static const std::string noPath{"<none>"};
        out << "n0 = " << options.n0
            << ", n1 = " << options.n1
            << ", f0 = " << options.f0
            << ", f1 = " << options.f1
            << ", output: `" << (options.output_path ? options.output_path->string()
            : noPath) << "'";
        return out;
    }
    

}
