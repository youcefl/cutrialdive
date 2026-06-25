/*
* MIT License
* Created on 2026.06.19
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#include "logger.hpp"

#include <iostream>

namespace cutrialdive {

    logger & log()
    {
        static logger log{verbosity::normal, std::cout};
        return log;
    }

}
