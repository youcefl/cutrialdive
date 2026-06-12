/*
* Creation date: 2026.04.29
* Created by Youcef Lemsafer
*/
#pragma once

#include <string>
#include <string_view>
#include <stdexcept>
#include <optional>
#include <charconv>

#include "num_seq_id.hpp"
#include "num_seq_spec.hpp"
#include "builtin_number_sequences.hpp"


namespace cutrialdive {
    template <typename IndexT, typename Func>
    auto dispatch_num_seq(num_seq_id numSeqId, Func&& func);    
}

// Impl
namespace cutrialdive {
    namespace details {
        ///
        template <typename NumberSequenceT, typename IndexT, typename Func, typename... Args>
        auto dispatch_num_seq_impl(Func&& func, Args&&...args)
        {
            return func.template operator()<NumberSequenceT>(std::forward<Args>(args)...);
        }

        template <typename IntT>
        std::optional<IntT> parse_int(std::string_view str)
        {
            IntT result{};
            auto [ptr, ec] = std::from_chars(str.data(), str.data() + str.size(), result);
            if(ec == std::errc()) {
                return {result};
            }
            return {};
        }

        template <typename IndexT, typename Func>
        auto dispatch_smarandache(uint64_t base, Func&& func)
        {
            switch(base) {
            case 2: return dispatch_num_seq_impl<smarandache<2>, IndexT>(func);
            case 3: return dispatch_num_seq_impl<smarandache<3>, IndexT>(func);
            case 4: return dispatch_num_seq_impl<smarandache<4>, IndexT>(func);
            case 5: return dispatch_num_seq_impl<smarandache<5>, IndexT>(func);
            case 6: return dispatch_num_seq_impl<smarandache<6>, IndexT>(func);
            case 7: return dispatch_num_seq_impl<smarandache<7>, IndexT>(func);
            case 8: return dispatch_num_seq_impl<smarandache<8>, IndexT>(func);
            case 9: return dispatch_num_seq_impl<smarandache<9>, IndexT>(func);
            case 10: return dispatch_num_seq_impl<smarandache<10>, IndexT>(func);
            case 11: return dispatch_num_seq_impl<smarandache<11>, IndexT>(func);
            case 12: return dispatch_num_seq_impl<smarandache<12>, IndexT>(func);
            case 13: return dispatch_num_seq_impl<smarandache<13>, IndexT>(func);
            case 14: return dispatch_num_seq_impl<smarandache<14>, IndexT>(func);
            case 15: return dispatch_num_seq_impl<smarandache<15>, IndexT>(func);
            case 16: return dispatch_num_seq_impl<smarandache<16>, IndexT>(func);
            case 17: return dispatch_num_seq_impl<smarandache<17>, IndexT>(func);
            case 18: return dispatch_num_seq_impl<smarandache<18>, IndexT>(func);
            case 19: return dispatch_num_seq_impl<smarandache<19>, IndexT>(func);
            case 20: return dispatch_num_seq_impl<smarandache<20>, IndexT>(func);
            case 21: return dispatch_num_seq_impl<smarandache<21>, IndexT>(func);
            case 22: return dispatch_num_seq_impl<smarandache<22>, IndexT>(func);
            case 23: return dispatch_num_seq_impl<smarandache<23>, IndexT>(func);
            case 24: return dispatch_num_seq_impl<smarandache<24>, IndexT>(func);
            case 25: return dispatch_num_seq_impl<smarandache<25>, IndexT>(func);
            case 26: return dispatch_num_seq_impl<smarandache<26>, IndexT>(func);
            case 27: return dispatch_num_seq_impl<smarandache<27>, IndexT>(func);
            case 28: return dispatch_num_seq_impl<smarandache<28>, IndexT>(func);
            case 29: return dispatch_num_seq_impl<smarandache<29>, IndexT>(func);
            case 30: return dispatch_num_seq_impl<smarandache<30>, IndexT>(func);
            case 31: return dispatch_num_seq_impl<smarandache<31>, IndexT>(func);
            case 32: return dispatch_num_seq_impl<smarandache<32>, IndexT>(func);
            case 33: return dispatch_num_seq_impl<smarandache<33>, IndexT>(func);
            case 34: return dispatch_num_seq_impl<smarandache<34>, IndexT>(func);
            case 35: return dispatch_num_seq_impl<smarandache<35>, IndexT>(func);
            case 36: return dispatch_num_seq_impl<smarandache<36>, IndexT>(func);
            case 37: return dispatch_num_seq_impl<smarandache<37>, IndexT>(func);
            case 38: return dispatch_num_seq_impl<smarandache<38>, IndexT>(func);
            case 39: return dispatch_num_seq_impl<smarandache<39>, IndexT>(func);
            case 40: return dispatch_num_seq_impl<smarandache<40>, IndexT>(func);
            case 41: return dispatch_num_seq_impl<smarandache<41>, IndexT>(func);
            case 42: return dispatch_num_seq_impl<smarandache<42>, IndexT>(func);
            case 43: return dispatch_num_seq_impl<smarandache<43>, IndexT>(func);
            case 44: return dispatch_num_seq_impl<smarandache<44>, IndexT>(func);
            case 45: return dispatch_num_seq_impl<smarandache<45>, IndexT>(func);
            case 46: return dispatch_num_seq_impl<smarandache<46>, IndexT>(func);
            case 47: return dispatch_num_seq_impl<smarandache<47>, IndexT>(func);
            case 48: return dispatch_num_seq_impl<smarandache<48>, IndexT>(func);
            case 49: return dispatch_num_seq_impl<smarandache<49>, IndexT>(func);
            case 50: return dispatch_num_seq_impl<smarandache<50>, IndexT>(func);
            case 51: return dispatch_num_seq_impl<smarandache<51>, IndexT>(func);
            case 52: return dispatch_num_seq_impl<smarandache<52>, IndexT>(func);
            case 53: return dispatch_num_seq_impl<smarandache<53>, IndexT>(func);
            case 54: return dispatch_num_seq_impl<smarandache<54>, IndexT>(func);
            case 55: return dispatch_num_seq_impl<smarandache<55>, IndexT>(func);
            case 56: return dispatch_num_seq_impl<smarandache<56>, IndexT>(func);
            case 57: return dispatch_num_seq_impl<smarandache<57>, IndexT>(func);
            case 58: return dispatch_num_seq_impl<smarandache<58>, IndexT>(func);
            case 59: return dispatch_num_seq_impl<smarandache<59>, IndexT>(func);
            case 60: return dispatch_num_seq_impl<smarandache<60>, IndexT>(func);
            case 61: return dispatch_num_seq_impl<smarandache<61>, IndexT>(func);
            case 62: return dispatch_num_seq_impl<smarandache<62>, IndexT>(func);
            case 63: return dispatch_num_seq_impl<smarandache<63>, IndexT>(func);
            case 64: return dispatch_num_seq_impl<smarandache<64>, IndexT>(func);
            }
            throw std::runtime_error{"Base `" + std::to_string(base) + "' is not supported for Smarandache numbers"};
        }
    } // namespace details


    template <typename IndexT, typename Func>
    auto dispatch_num_seq(num_seq_spec numSeqSpec, Func&& func)
    {
        using details::dispatch_num_seq_impl;
        auto numSeqId = numSeqSpec.seq_id;
        switch (numSeqId) {
            case num_seq_id::mersenne:
                return dispatch_num_seq_impl<mersenne, IndexT>(func);
            case num_seq_id::proth: {
                auto k = details::parse_int<uint64_t>(numSeqSpec.seq_params);
                if(!k) {
                    throw std::runtime_error{"Missing or invalid k value `" + numSeqSpec.seq_params + "' in Proth number spec"};
                }
                return dispatch_num_seq_impl<proth, IndexT>(func, *k);
            }
            case num_seq_id::smarandache: {
                auto base = numSeqSpec.seq_params.empty()
                            ? std::optional<uint64_t>{10}
                            : details::parse_int<uint64_t>(numSeqSpec.seq_params);
                if(!base) {
                    throw std::runtime_error{"Invalid base `" + numSeqSpec.seq_params + "' for Smarandache numbers"};
                }
                return details::dispatch_smarandache<IndexT>(*base, func);
            }
            default: throw details::unexpected_num_seq_id_exception(numSeqId);
        }
    }

}
