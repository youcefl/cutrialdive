/*
* Created on 2026.06.10
* Copyright (c) Youcef Lemsafer
*/
#pragma once

#include <iostream>

namespace cutrialdive {

    struct trial_factoring_options;
    struct tf_runtime_options;
    template <typename T, typename U>
    class factoring_results;
    class checkpoint_manager;
    struct engine_state;


    struct trial_factoring_context
    {
        trial_factoring_options const & options;
        tf_runtime_options const & runtime_options;
        factoring_results<uint64_t, uint32_t> & results;
        std::ostream & output_stream;
        checkpoint_manager * checkpoint{};
        engine_state * resume_state{};
    };
}
