/*
* MIT License
* Created on 2026.07.13
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#include "fibonacci.hpp"
#include "modular_arithmetic_detail.hpp"


namespace cutrialdive {

    HgInt fibonacci::value(index_type n)
    {
        HgInt nextValue{};
        return value(n, nextValue);
    }

    std::ostream& fibonacci::print_value(index_type n, std::ostream & out)
    {
        return out << value(n);
    }

    char const* fibonacci::short_name()
    {
        return "I";
    }

    HgInt fibonacci::value(index_type n, HgInt & nextValue)
    {
        HgInt f0{};
        nextValue = HgInt{1};
        HgInt & f1 = nextValue;
        if(!n) {
            return f0;
        }
        for(index_type m = index_type{1} << (bit_length(n) - 1); m; m >>= 1) {
            if(m & n) {
                auto t0 = (2*f1 - f0)*f0;
                f0 = f1*f1 + f0*f0;
                f1 = f0 + t0;
            } else {
                auto t0 = f0;
                f0 = (2*f1 - f0)*f0;
                f1 = f1*f1 + t0*t0;
            }
        }
        return f0;
    }

}
