/*
* Creation date: 2026.04.29
* Created by Youcef Lemsafer
*/
#include "mode.hpp"

#include <unordered_map>

namespace cutrialdive {

    mode_flag mode_from_string(std::string const & modeAsStr)
    {
        static const std::unordered_map<std::string, mode_flag> flags_map {
            {"mersenne", mode_flag::mersenne},
            {"smarandache", mode_flag::smarandache},
        }
        auto it = flags_map.find(modeAsStr);
        if(it == flags_map.end()) {
            throw std::runtime_error{std::string{"Invalid mode `"} + modeAsStr + "'"};
        }
        return it-> second;
    }

}
