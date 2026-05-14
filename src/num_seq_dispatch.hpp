/*
* Creation date: 2026.04.29
* Created by Youcef Lemsafer
*/
#pragma once

#include <string>
#include <stdexcept>
#include "num_seq_id.hpp"
#include "number_sequence.hpp"

namespace cutrialdive {


    template <typename IndexT, typename Func>
    auto dispatch_num_seq(num_seq_id numSeqId, Func&& func);    
}

namespace cutrialdive {
    namespace details {

        std::runtime_error unexpected_num_seq_id_exception(num_seq_id numSeqId);

        template <num_seq_id numSeqId, typename IndexT, typename Func>
        auto dispatch_num_seq_impl(Func&& func)
        {
            return func.template operator()<number_sequence<numSeqId>>();
        }
    } // namespace details


    template <typename IndexT, typename Func>
    auto dispatch_num_seq(num_seq_id numSeqId, Func&& func)
    {
        using details::dispatch_num_seq_impl;
        switch (numSeqId) {
            case num_seq_id::mersenne:    return dispatch_num_seq_impl<num_seq_id::mersenne,    IndexT>(func);
            case num_seq_id::smarandache: return dispatch_num_seq_impl<num_seq_id::smarandache, IndexT>(func);
            default: throw details::unexpected_num_seq_id_exception(numSeqId);
        }
    }

}
