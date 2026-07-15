/*
* MIT License
* Created on 2026.07.13
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#include "lucas.hpp"

namespace cutrialdive {

    HgInt lucas::value(index_type n)
    {
        HgInt fn1; // F_{n+1}
        // F_n
        HgInt fn = fibonaccci_.value(n, fn1);
        // L_n = 2F_{n+1} - F_n
        return 2*fn1 - fn;
    }

    std::ostream& lucas::print_value(index_type n, std::ostream & out)
    {
        return out << value(n);
    }

    char const* lucas::short_name()
    {
        return "L";
    }

}
