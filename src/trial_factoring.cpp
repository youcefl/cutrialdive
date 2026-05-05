/*
* Created on 2026.05.04
* Copyright (c) Youcef Lemsafer
*/
#include "trial_factoring.hpp"

#include <iostream>
#include <fstream>
#include <omp.h>

#include "factors_buffer.hpp"
#include "factoring_results.hpp"
#include "hgint.hpp"
#include "siever.hpp"
#include "timer.hpp"
#include "modular_arithmetic_detail.hpp"
#include "number_sequence.hpp"
#include "builtin_number_sequences.hpp"



namespace {

    using namespace cutrialdive;

    template <typename PrimeT, typename ExponentT>
    void read_previous_results(uint64_t n0, uint64_t n1, factoring_results<PrimeT, ExponentT> & previous_results)
    {
        if(n1 < n0) {
            throw std::logic_error{"Invalid range specification, n0 must be less than or equal to n1."};
        }
        if(n0 != previous_results.n0() || n1 - n0 != previous_results.size()) {
            std::ostringstream ostr;
            ostr << "Mismatch between factoring results: previous(n0=" 
                    << previous_results.n0() << ", n1=" << previous_results.n0() + previous_results.size()
                    << "), current(n0=" << n0 << ", n1=" << n1 << ").";
            throw std::runtime_error{ostr.str()};
        }

        // @todo: Implement update of previous results
        throw std::runtime_error{"Update of previous results not implemented yet"};
    }

}

namespace cutrialdive {

    namespace {
        // Will we ever have a trial factoring reveal more than 128 factors in normal usage?
        uint32_t MAX_FACTORS_PER_NUMBER = 128;
    }

    template <typename NumberSequenceT>
        requires NumberSequence<NumberSequenceT>
    void do_trial_factor(
        trial_factoring_options const & opts,
        factoring_results<uint64_t, uint32_t> & results
        )
    {
        constexpr auto max_num_of_primes = size_t{1} << 22;
        constexpr auto segment_len = uint64_t{1} << 26;

        uint64_t n0 = opts.n0, n1 = opts.n1, f0 = opts.f0, f1 = opts.f1;

        std::vector<uint64_t> primes, inv_primes;
        primes.reserve(max_num_of_primes);
        inv_primes.reserve(max_num_of_primes);
        // residue of S(n0) modulo each prime
        std::vector<uint64_t> n0_residues;
        factors_buffer<uint64_t, uint32_t> factors_buf{n0, n1 - n0, MAX_FACTORS_PER_NUMBER};

        auto const sm_n0 = NumberSequenceT::value(n0);

        auto init_time = std::chrono::duration<double>{},
             propagation_time = init_time;

        for(auto fx = f0, fy = (std::min)(fx + segment_len, f1);
            fy <= f1;
            fx = fy, fy = (std::max)(fx + segment_len, f1)
        ) {
            primes.clear();
            inv_primes.clear();
            sieve(fx, fy, primes);
            if(n0_residues.size() < primes.size()) {
                n0_residues.reserve(primes.size());
                n0_residues.assign(primes.size(), uint64_t{});
            } else if(n0_residues.size() > primes.size()) {
                n0_residues.resize(primes.size());
            }
            inv_primes.assign(primes.size(), 0);

            // Compute S(n0) mod p and, if p != 2, the inverse of p mod 2^64.
            auto startTime = std::chrono::high_resolution_clock::now();
            #pragma omp parallel
            {
                HgInt local_sm_n0{sm_n0};
                size_t const i_max = primes.size();
                #pragma omp for
                for(size_t i = 0; i < i_max; ++i) {
                    auto p = primes[i];
                    n0_residues[i] = local_sm_n0.mod(p);
                    if(!n0_residues[i]) {
                        factors_buf.push_factor(0, p);
                    }
                    // We still push a value when p is 2 to avoid having the indice of
                    // p and mu(p) differ by one.
                    inv_primes[i] = p == 2 ? 0 : mu(p);
                }
            }
            auto now = std::chrono::high_resolution_clock::now();
            init_time += std::chrono::duration<double>(now - startTime);
            startTime = now;
            auto n = n0;
            for(auto i = 1; i < n1 - n0; ++i, ++n) {
                // 2 is so special that it dictates a special treatment.
                // Sieve results are always ordered so 2 can only be at the front.
                bool have_to_skip_front = false;
                if(!primes.empty() && primes.front() == 2) {
                    n0_residues[0] = NumberSequenceT::next_value_mod_2(n0_residues[0], n);
                    if(!n0_residues[0]) {
                        factors_buf.push_factor(i, 2);
                    }
                    have_to_skip_front = true;
                }
                size_t j0 = have_to_skip_front ? 1 : 0, j_max = primes.size();
                #pragma omp parallel for
                for(size_t j = j0;
                    j < j_max;
                    ++j) {
                    n0_residues[j] = NumberSequenceT::next_value_mod_mu(n0_residues[j], n, primes[j], inv_primes[j]);
#if 0
                    auto a = __uint128_t(n0_residues[j]) * 100000 + n;
                    auto q = (a * inv_primes[j]) >> 64;
                    n0_residues[j] = a - q * primes[j];
                    if(n0_residues[j] >= primes[j]) {
                        n0_residues[j] -= primes[j];
                    }
#endif
                    if(!n0_residues[j]) {
                        factors_buf.push_factor(i, primes[j]);
                    }
                }
            }
            now = std::chrono::high_resolution_clock::now();
            propagation_time += std::chrono::duration<double>(now - startTime);
        }
        std::cout << "Initialization time: " << init_time.count() << "s" << std::endl;
        std::cout << "Propagation time: " << propagation_time.count() << "s" << std::endl;
        results = factors_buf.to_factoring_results();
    }

    template <typename NumberSequenceT>
        requires NumberSequence<NumberSequenceT>
    void trial_factor(
        trial_factoring_options const & opts,
        factoring_results<uint64_t, uint32_t> & previous_results
    )
    {
        if(opts.n1 <= opts.n0) {
            std::cout << "Range [" << opts.n0 << ", " << opts.n1 << "[ is empty, no number to process..." << std::endl;
            return;
        }
        if(opts.f1 <= opts.f0) {
            std::cout << "Range [" << opts.f0 << ", " << opts.f1 << "[ is empty, no prime numbers to divide by..." << std::endl;
            return;
        }
        if(previous_results.empty()) {
            previous_results.reserve(opts.n1 - opts.n0);
        } else {
            read_previous_results(opts.n0, opts.n1, previous_results);
        }
        do_trial_factor<NumberSequenceT>(opts, previous_results);
    }

    /// @brief Performs trial factoring according to given options
    /// @param opts trial factoring options
    void trial_factor(
        mode_flag modeFlag,
        trial_factoring_options const & opts
    )
    {
        std::unique_ptr<std::ofstream> output_ptr;
        if(opts.output_path) {
            auto path = opts.output_path.value();
            if(std::filesystem::exists(path)) {
                std::string msg{"Output: cannot overwrite existing file "};
                msg += path.string();
                throw std::runtime_error(msg);
            }
            output_ptr = std::make_unique<std::ofstream>(path);
        }
        factoring_results<uint64_t, uint32_t> results{opts.n0, opts.n1 - opts.n0};
        std::cout << "Trial factoring starts, results will be written to "
            << (output_ptr ? (std::string{"file `"} + opts.output_path.value().string() + "'")
                           : "the console")
            << std::endl;
        std::cout << "Maximum number of threads: " << omp_get_max_threads() << std::endl;
        {
            timer tfTimer{"Trial factoring took ", std::cout};
            
            dispatch_mode<decltype(opts.n0)>(modeFlag, [&]<typename Seq>() {
                trial_factor<Seq>(opts, results);
            });
        }
        (output_ptr ? *output_ptr.get() : std::cout) << results;
    }


}
