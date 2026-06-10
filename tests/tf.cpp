/*
* Created on 2026.06.04
* Copyright (c) Youcef Lemsafer
*/
#include "tf.hpp"

#include <sstream>

#include "peano.hpp"
#include "trial_factoring_impl.hpp"

namespace cutrialdive::tests {

    template <typename PeanoT>
    factoring_results<uint64_t, uint32_t> tf_peano(trial_factoring_options const & tfOptions, std::ostream & out)
    {
        factoring_results<uint64_t, uint32_t> results{tfOptions.n0, tfOptions.n1 - tfOptions.n0};
        trial_factor<PeanoT>(tfOptions,tf_runtime_options::default_options(), results, out);
        return results;
    }

    // peano {value_mod, next_value_mod}
    factoring_results<uint64_t, uint32_t> basic_tf_peano(trial_factoring_options const & tfOptions, std::ostream & out)
    {
        return tf_peano<peano<function_to_define::value_mod,
                        function_to_define::next_value_mod>>(tfOptions, out);
    }

    // peano {value, value_mod_2, value_mod_mu64, next_value_mod}
    factoring_results<uint64_t, uint32_t> value_mod_2_mu64_tf_peano(trial_factoring_options const & tfOptions, std::ostream & out)
    {
        return tf_peano<peano<function_to_define::value,
                        function_to_define::value_mod_2,
                        function_to_define::value_mod_mu64,
                        function_to_define::next_value_mod>>(tfOptions, out);
    }

    // peano {value, value_mod_2, value_mod_mu128, next_value_mod_2, next_value_mod_mu128}
    factoring_results<uint64_t, uint32_t> value_mod_2_mu128_tf_peano(trial_factoring_options const & tfOptions, std::ostream & out)
    {
        return tf_peano<peano<function_to_define::value,
                        function_to_define::value_mod_2,
                        function_to_define::value_mod_mu128,
                        function_to_define::next_value_mod_2,
                        function_to_define::next_value_mod_mu128>>(tfOptions, out);
    }

}
