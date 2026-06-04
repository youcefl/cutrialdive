/*
* Created on 2026.06.04
* Copyright (c) Youcef Lemsafer
*/
#pragma once

#include <ostream>
#include "trial_factoring_options.hpp"
#include "factoring_results.hpp"

namespace cutrialdive::tests {

    // peano {value_mod, next_value_mod}
    factoring_results<uint64_t, uint32_t> basic_tf_peano(trial_factoring_options const & tfOptions, std::ostream & out);

    // peano {value, value_mod_2, value_mod_mu64, next_value_mod}
    factoring_results<uint64_t, uint32_t> value_mod_2_mu64_tf_peano(trial_factoring_options const & tfOptions, std::ostream & out);

    // peano {value, value_mod_2, value_mod_mu128, next_value_mod_2, next_value_mod_mu128}
    factoring_results<uint64_t, uint32_t> value_mod_2_mu128_tf_peano(trial_factoring_options const & tfOptions, std::ostream & out);
}
