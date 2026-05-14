/*
* Creation date: 2026.04.29
* Created by Youcef Lemsafer
*/
#include "mode.hpp"

#include <sstream>
#include <unordered_map>
#include <stdexcept>

namespace cutrialdive {
    namespace details {

        std::runtime_error unexpected_mode_exception(mode_flag modeFlag)
        {
            std::ostringstream ostr;
            ostr << "Unexpected mode " << static_cast<int>(modeFlag);
            return std::runtime_error{ostr.str()};
        }

    } // namespace details

    mode_flag mode_from_string(std::string const & modeAsStr)
    {
        static const std::unordered_map<std::string, mode_flag> flags_map {
            {"mersenne", mode_flag::mersenne},
            {"smarandache", mode_flag::smarandache},
        };
        auto it = flags_map.find(modeAsStr);
        if(it == flags_map.end()) {
            throw std::runtime_error{std::string{"Invalid mode `"} + modeAsStr + "'"};
        }
        return it-> second;
    }


}
