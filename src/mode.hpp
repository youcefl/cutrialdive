/*
* Creation date: 2026.04.29
* Created by Youcef Lemsafer
*/
#pragma once

#include <string>
#include <stdexcept>
#include "mode_flag.hpp"
#include "number_sequence.hpp"

namespace cutrialdive {

    /// @brief Returns the mode as a cutrialdive::mode value
    mode_flag mode_from_string(std::string const & modeAsStr);

    template <typename IndexT, typename Func>
    auto dispatch_mode(mode_flag modeFlag, Func&& func);    
}

namespace cutrialdive {
    namespace details {

        std::runtime_error unexpected_mode_exception(mode_flag modeFlag);

        template <mode_flag modeFlag, typename IndexT, typename Func>
        auto dispatch_mode_impl(Func&& func)
        {
            return func.template operator()<number_sequence<modeFlag>>();
        }
    } // namespace details


    template <typename IndexT, typename Func>
    auto dispatch_mode(mode_flag modeFlag, Func&& func)
    {
        using details::dispatch_mode_impl;
        switch (modeFlag) {
            case mode_flag::mersenne:    return dispatch_mode_impl<mode_flag::mersenne,    IndexT>(func);
            case mode_flag::smarandache: return dispatch_mode_impl<mode_flag::smarandache, IndexT>(func);
            default: throw details::unexpected_mode_exception(modeFlag);
        }
    }

}
