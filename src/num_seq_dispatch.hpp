/*
* Creation date: 2026.04.29
* Created by Youcef Lemsafer
*/
#pragma once

#include <string>
#include <stdexcept>
#include "num_seq_id.hpp"
#include "num_seq_spec.hpp"
#include "number_sequence.hpp"

namespace cutrialdive {
    template <typename IndexT, typename Func>
    auto dispatch_num_seq(num_seq_id numSeqId, Func&& func);    
}

// Impl
namespace cutrialdive {
    namespace details {
        ///
        template <typename NumberSequenceT, typename IndexT, typename Func>
        auto dispatch_num_seq_impl(Func&& func)
        {
            return func.template operator()<NumberSequenceT>();
        }
    } // namespace details


    template <typename IndexT, typename Func>
    auto dispatch_num_seq(num_seq_spec numSeqSpec, Func&& func)
    {
        using details::dispatch_num_seq_impl;
        auto numSeqId = numSeqSpec.num_seq_id_value;
        switch (numSeqId) {
            case num_seq_id::mersenne:
                return dispatch_num_seq_impl<number_sequence<num_seq_id::mersenne>, IndexT>(func);
            case num_seq_id::smarandache: {
                if(numSeqSpec.num_seq_params.empty() || numSeqSpec.num_seq_params == "10") {
                    return dispatch_num_seq_impl<number_sequence<num_seq_id::smarandache>, IndexT>(func);
                }
                if(numSeqSpec.num_seq_params == "2") {
                    return dispatch_num_seq_impl<number_sequence<num_seq_id::smarandache, 2>, IndexT>(func);
                }
                if(numSeqSpec.num_seq_params == "3") {
                    return dispatch_num_seq_impl<number_sequence<num_seq_id::smarandache, 3>, IndexT>(func);
                }
                throw std::runtime_error{"Invalid base `" + numSeqSpec.num_seq_params + "' for Smarandache numbers"};
            }
            default: throw details::unexpected_num_seq_id_exception(numSeqId);
        }
    }

}
