/*
* Youcef Lemsafer, 2026.04.21
* Definition of the sieving back-end
*/
#pragma once

#include <vector>
#include <cstdint>
#ifdef CUTRIALDIVE_HAS_PRIMESIEVE
#include <primesieve.hpp>
#endif
#ifdef CUTRIALDIVE_HAS_LFP
#include <lfp.hpp>
#include <omp.h>
#endif

#if !defined CUTRIALDIVE_HAS_PRIMESIEVE && !defined CUTRIALDIVE_HAS_LFP
#error "One of CUTRIALDIVE_HAS_PRIMESIEVE and CUTRIALDIVE_HAS_LFP must be defined!"
#endif

namespace cutrialdive {

    enum sieve_lib
    {
        primesieve,
        lfp
    };

    template <uint32_t SieveLib>
    void sievex(uint64_t start, uint64_t end, std::vector<uint64_t> & primes)
    {
        static_assert((SieveLib == sieve_lib::primesieve)
                     || (SieveLib == sieve_lib::lfp));

        if constexpr (SieveLib == sieve_lib::primesieve) {
#ifdef CUTRIALDIVE_HAS_PRIMESIEVE
            primesieve::generate_primes(start, end ? end - 1 : end, &primes);
#else
            static_assert(SieveLib != sieve_lib::primesieve, "Sieving library primesieve is not available");
#endif
        }
        else if constexpr(SieveLib == sieve_lib::lfp) {
#ifdef CUTRIALDIVE_HAS_LFP
            auto sieveRes = lfp::sieve<uint64_t>(start, end, lfp::threads{uint32_t(omp_get_max_threads())});
            auto rng = sieveRes.range();
            primes.assign(rng.begin(), rng.end());
#else
            static_assert(SieveLib != sieve_lib::lfp, "Sieving library lfp is not available");
#endif
        }
    }

#ifdef CUTRIALDIVE_HAS_PRIMESIEVE
    inline void sieve(uint64_t start, uint64_t end, std::vector<uint64_t> & primes)
    {
        sievex<sieve_lib::primesieve>(start, end, primes);
    }
#elif defined CUTRIALDIVE_HAS_LFP 
    inline void sieve(uint64_t start, uint64_t end, std::vector<uint64_t> & primes)
    {
        sievex<sieve_lib::lfp>(start, end, primes);
    }
#endif
}
