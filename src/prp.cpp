/*
 * Youcef Lemsafer
 * 2023.01.19
 */
#include "hgint.hpp"
#include "gw_utility.h"

bool isPRP(HgInt number)
{
    return gw_prp(number.get()) != 0;
}

bool is_prp(HgInt number)
{
    return isPRP(number);
}

