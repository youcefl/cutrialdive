/*
* Creation date: 2026.05.05
* Created by Youcef Lemsafer
*/
#pragma once

#include <cstdint>
#include <vector>
#include <ostream>

#include "mode_flag.hpp"
#include "factor.hpp"


namespace cutrialdive {

    void run_prp_test(
        mode_flag modeFlag,
        uint64_t n,
        std::vector<factor<uint64_t, uint32_t>> & factors,
        bool haveToBoostFactors,
        std::ostream & out
    );

}

