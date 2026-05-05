/*
* Creation date: 2026.05.05
* Created by Youcef Lemsafer
*/
#pragma once

namespace cutrialdive {

    /// Available modes
    enum class mode_flag {
        mersenne,
        smarandache
    };

    template <typename StrT>
    StrT to_string(mode_flag modeFlag);
}

namespace cutrialdive {

    template <typename StrT>
    StrT to_string(mode_flag modeFlag)
    {
        switch(modeFlag) {
            case mode_flag::mersenne: return {"mersenne"};
            case mode_flag::smarandache: return {"smarandache"};
            default: "unknown";
        }
    }

}
