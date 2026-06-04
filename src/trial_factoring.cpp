/*
* Created on 2026.05.04
* Copyright (c) Youcef Lemsafer
*/
#include "trial_factoring.hpp"
#include "trial_factoring_impl.hpp"

namespace cutrialdive {

    /// @brief Performs trial factoring according to given options
    /// @param opts trial factoring options
    factoring_results<uint64_t, uint32_t>
    trial_factor(
        num_seq_spec numSeqSpec,
        trial_factoring_options const & opts,
        std::ostream & out
    )
    {
        factoring_results<uint64_t, uint32_t> results{opts.n0, opts.n1 - opts.n0};
        dispatch_num_seq<decltype(opts.n0)>(numSeqSpec, [&]<typename Seq>() {
            trial_factor<Seq>(opts, results, out);
        });
        return results;
    }
}
