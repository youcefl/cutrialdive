/*
* Created on 2026.06.19
* Copyright (c) Youcef Lemsafer
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
