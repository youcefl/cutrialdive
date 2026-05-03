/* Created on 2026.05.03
* Copyright (c) Youcef Lemsafer
*/


// Below is an attempt at making a tiled residue computation: load tile of N
// from global memory into shared memory, compute residue move on to next tile
// Did not really improve performance (less than 2%) and forces us to
// preseve residues across tiles.
// Abandoned for the moment...

template <uint8_t BatchSize, typename PrimeT>
__global__
void trial_div2(uint64_t * number, uint64_t numberLen,
    typename device_prime_data<PrimeT>::data primeData,
    int* deviceFactorsCount,
    uint64_t * deviceFactors,
    uint32_t deviceFactorsSize,
    uint64_t * residues,
    uint64_t residuesSize
)
{
    static_assert((BatchSize > 0) && (BatchSize <= 3), "The only batch sizes supported are 1, 2 and 3");

//     if(threadIdx.x == 0 && blockIdx.x == 0) {
//         printf(R"-(Kernel called
//     numberLen = %llu
//     primes count: %llu
//     device factors count: %d
//     device factors max: %u
//     residues size: %llu
// )-", numberLen, primeData.primes_count, *deviceFactorsCount, deviceFactorsSize, residuesSize);
//     }


    constexpr unsigned int maxTileSize = 4096;
    __shared__ uint64_t numberTile[maxTileSize];

    // Copy tiles into shared memory starting from most significant limbs
    // going backwards till less significant limbs
    for(auto numIdxEnd = numberLen, numIdxStart = (numIdxEnd >= maxTileSize) ? numIdxEnd - maxTileSize : 0;
        numIdxEnd > 0;
        numIdxEnd = numIdxStart,
        numIdxStart = (numIdxStart >= maxTileSize) ? numIdxStart - maxTileSize : 0
    ) {
        auto numBaseIdx = numIdxStart + threadIdx.x;
        auto numStride = blockDim.x;
        // if(threadIdx.x == 0 && blockIdx.x == 0) {
        //     printf("Copying bytes of N [%llu, %llu[\n", numIdxStart, numIdxEnd);
        // }

        for(auto i = numBaseIdx; i < numIdxEnd; i += numStride) {
            if(i - numIdxStart >= maxTileSize) {
                break;
            }
            numberTile[i - numIdxStart] = number[i];
        }

        // Wait for all threads to be done with the tile copy
        __syncthreads();

        const auto tileSize = numIdxEnd - numIdxStart;
        // if(threadIdx.x == 0 && blockIdx.x == 0) {
        //     printf("Tile size: %llu\n", tileSize);
        // }
        const bool isLastTile = numIdxStart == 0;

        // Perform divisions on current tile
        int index = (threadIdx.x + blockIdx.x * blockDim.x) * BatchSize;
        int stride = blockDim.x * gridDim.x * BatchSize;

        auto p = primeData.primes;
        auto plen = primeData.primes_count;
        for(size_t i = index; i < plen; i += stride) {
            if constexpr (BatchSize == 1) {
                residues[i] = modnby1(numberTile, tileSize, p[i], primeData.reciprocals[i],
                                        initial_residue<PrimeT>{residues[i]});
                if(isLastTile && !residues[i]) {
                    pushFactor(p[i], deviceFactorsCount, deviceFactors, deviceFactorsSize);
                }
            } else {
                uint64_t pp = p[i];
                #pragma unroll
                for(auto j = 1; j < BatchSize; ++j) {
                    if (i + j >= plen) {
                        break;
                    }
                    pp *= p[i + j];
                }
                residues[i] = modnby1(numberTile, tileSize, pp, initial_residue<PrimeT>{residues[i]});
                if(isLastTile) {
                    #pragma unroll
                    for (auto j = 0; j < BatchSize; ++j) {
                        if (i + j >= plen) {
                            break;
                        }
                        residues[i + j] = residues[i] % p[i + j];
                        if(isLastTile && !residues[i + j]) {
                            pushFactor(p[i + j], deviceFactorsCount, deviceFactors, deviceFactorsSize);
                        }
                    }
                }
            }
        }

        // Wait for all threads to be done with the division
        __syncthreads();
    }
}
