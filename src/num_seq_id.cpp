/*
* Creation date: 2026.04.29
* Created by Youcef Lemsafer
*/
#include "num_seq_id.hpp"

#include <sstream>
#include <unordered_map>
#include <stdexcept>

namespace cutrialdive {
    namespace details {

        std::runtime_error unexpected_num_seq_id_exception(num_seq_id numSeqId)
        {
            std::ostringstream ostr;
            ostr << "Unexpected number sequence id: " << static_cast<int>(numSeqId);
            return std::runtime_error{ostr.str()};
        }

    } // namespace details

    num_seq_id num_seq_id_from_string(std::string const & numSeqIdAsStr)
    {
        static const std::unordered_map<std::string, num_seq_id> flags_map {
            {"mersenne", num_seq_id::mersenne},
            {"smarandache", num_seq_id::smarandache},
        };
        auto it = flags_map.find(numSeqIdAsStr);
        if(it == flags_map.end()) {
            throw std::runtime_error{std::string{"Invalid number sequence id `"} + numSeqIdAsStr + "'"};
        }
        return it-> second;
    }


}
