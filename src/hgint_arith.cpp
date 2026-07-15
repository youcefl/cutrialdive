/*
* MIT License
* Created on 2026.07.15
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#include "hgint_arith.hpp"

namespace cutrialdive {

    hgint_arith::value_type hgint_arith::mul_by_2(value_type a) const
    {
        return a << 1;
    }
    hgint_arith::value_type hgint_arith::add(value_type a, value_type b) const
    {
        return a + b;
    }
    hgint_arith::value_type hgint_arith::sub(value_type a, value_type b) const
    {
        return a - b;
    }
    hgint_arith::value_type hgint_arith::mul(value_type a, value_type b) const
    {
        return a * b;
    }
    hgint_arith::value_type hgint_arith::sqr(value_type a) const
    {
        return a * a;
    }

}
