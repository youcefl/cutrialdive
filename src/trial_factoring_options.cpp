/*
* MIT License
* Created on 2026.05.05
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
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
                                                     : noPath)
            << "'";
        return out;
    }

    tf_runtime_options tf_runtime_options::default_options()
    {
        return tf_runtime_options{
            .common_options = runtime_options::default_options(),
            .checkpoint_period = std::chrono::seconds{360},
            .progress_period = std::chrono::milliseconds{1000},
            .is_progress_enabled = true,
            .segment_length = uint64_t{1} << 26,
#ifdef CUTRIALDIVE_ENABLE_GPU
            .device_id = 0
#endif
        };
    }

}
