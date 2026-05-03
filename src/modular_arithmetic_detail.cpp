/*
* Created on 2026.04.30
* Copyright Youcef Lemsafer
*/
#include "modular_arithmetic_detail.hpp"

#include <iostream>
#include <omp.h>

namespace cutrialdive {

    void compute_reciprocals(std::vector<uint64_t> const & primes, std::vector<uint64_t> & reciprocalsVec)
    {
        if(primes.empty()) {
            reciprocalsVec.clear();
            return;
        }
        uint64_t const primes_size = primes.size();
        if(primes_size > reciprocalsVec.size()) {
            reciprocalsVec.resize(primes_size);
        }
        reciprocalsVec.assign(primes_size, 0);

        auto * reciprocals = &reciprocalsVec[0];

        constexpr auto blockSize = 4;


        if(primes_size >= blockSize) {
            #pragma omp parallel for schedule(static)
            for(uint64_t i = 0; i < primes_size - (blockSize - 1); i += blockSize) {
                reciprocals[i] = reciprocal(normalize(primes[i]));
                reciprocals[i] = reciprocal(normalize(primes[i+1]));
                reciprocals[i] = reciprocal(normalize(primes[i+2]));
                reciprocals[i] = reciprocal(normalize(primes[i+3]));
            }
        }

        for(uint64_t i = primes_size - (primes_size % 4); i < primes_size; ++i) {
            reciprocals[i] = reciprocal(normalize(primes[i]));
        }
    }


}