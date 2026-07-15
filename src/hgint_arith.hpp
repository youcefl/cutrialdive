/*
* MIT License
* Created on 2026.07.15
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#pragma once


#include <cstdint>
#include "hgint.hpp"


namespace cutrialdive {

    struct hgint_arith
    {
        using value_type = HgInt;
        value_type mul_by_2(value_type a) const;
        value_type add(value_type a, value_type b) const;
        value_type sub(value_type a, value_type b) const;
        value_type mul(value_type a, value_type b) const;
        value_type sqr(value_type a) const;
    };

}
