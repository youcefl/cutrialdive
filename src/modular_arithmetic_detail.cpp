/*
* Created on 2026.04.30
* Copyright Youcef Lemsafer
*/
#include "modular_arithmetic_detail.hpp"

#include <iostream>
#include <numeric>
#include <omp.h>

namespace cutrialdive {

    void compute_reciprocals(std::vector<std::span<uint64_t>> const & primes, std::vector<uint64_t> & reciprocalsVec)
    {
        if(primes.empty()) {
            reciprocalsVec.clear();
            return;
        }
        uint64_t const primes_size = std::accumulate(std::begin(primes), std::end(primes),
            uint64_t{}, [](uint64_t x, auto const & y) {
                return x + y.size();
            });
        if(primes_size > reciprocalsVec.size()) {
            reciprocalsVec.resize(primes_size);
        }
        reciprocalsVec.assign(primes_size, 0);

        constexpr auto blockSize = 4;

        auto j = 0;
        for(auto chunk : primes) {
            auto * reciprocals = &reciprocalsVec[j];
            auto chunkData = chunk.data();

            auto const chunkSize = chunk.size();
            if(chunkSize >= blockSize) {
                #pragma omp parallel for schedule(static)
                for(uint64_t i = 0; i < chunkSize - (blockSize - 1); i += blockSize) {
                    reciprocals[i] = reciprocal(normalize(chunkData[i]));
                    reciprocals[i] = reciprocal(normalize(chunkData[i+1]));
                    reciprocals[i] = reciprocal(normalize(chunkData[i+2]));
                    reciprocals[i] = reciprocal(normalize(chunkData[i+3]));
                }
            }

            for(uint64_t i = chunkSize - (chunkSize % 4); i < chunkSize; ++i) {
                reciprocals[i] = reciprocal(normalize(chunkData[i]));
            }
            j += chunkSize;
        }
    }


}