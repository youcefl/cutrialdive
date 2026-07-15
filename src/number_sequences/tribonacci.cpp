/*
* MIT License
* Created on 2026.07.14
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#include "tribonacci.hpp"
#include "hgint_arith.hpp"

namespace cutrialdive {

    HgInt tribonacci::value(index_type n)
    {
        HgInt tn, tn1, tn2;
        tribonacci_at(hgint_arith{}, n, tn, tn1, tn2);
        return tn;
    }

    std::ostream& tribonacci::print_value(index_type n, std::ostream & out)
    {
        return out << value(n);
    }

    char const * tribonacci::short_name()
    {
        return "T";
    }

}
