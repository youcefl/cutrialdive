/*
* Created on 2026.06.04
* Copyright (c) Youcef Lemsafer
*/
#pragma once

#include <cstdint>
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h>

#ifdef CUTRIALDIVE_ENABLE_GPU
#include "trial_factoring.cuh"
#endif

#include "factors_buffer.hpp"
#include "factoring_results.hpp"
#include "hgint.hpp"
#include "siever.hpp"
#include "timer.hpp"
#include "modular_arithmetic_detail.hpp"
#include "number_sequence.hpp"
#include "number_sequence_helpers.hpp"
#include "num_seq_dispatch.hpp"
#include "barrett_reciprocals.hpp"
#include "builtin_number_sequences.hpp"
#include "progress.hpp"


namespace cutrialdive {
namespace detail {

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

} // namespace cutrialdive::detail


    template <typename NumberSequenceT>
        requires NumberSequence<NumberSequenceT>
    void host_trial_factor(
        trial_factoring_options const & opts,
        factoring_results<uint64_t, uint32_t> & results,
        std::ostream & out
        )
    {
        uint64_t n0 = opts.n0, n1 = opts.n1, f0 = opts.f0, f1 = opts.f1;

        factors_buffer<uint64_t, uint32_t> factors_buf{n0, n1 - n0, opts.max_factors_per_number};

        NumberSequenceT numSeq;

        auto progressHandler = opts.is_progress_enabled ? std::make_unique<progress>(f1, out)
                                                        : std::unique_ptr<progress>{};
        auto updateProgress = progressHandler
            ? std::function<void(uint64_t, uint64_t)>{[&progressHandler](auto segUpperBound, auto lastPrimeInSeg) {
                progressHandler->update(segUpperBound, lastPrimeInSeg); 
              }}
            : std::function<void(uint64_t, uint64_t)>{[](auto, auto) {}};

        // 2 is a special case, deal with it now so that only odd primes remain
        if(f0 <= 2 && f1 > 2) {
            f0 = 3;
            auto residue = value_mod_2(numSeq, n0);
            if(!residue) {
                factors_buf.push_factor(0, 2);
            }
            for(auto n = n0, nEnd = n1 - 1; n < nEnd; ++n) {
                residue = next_value_mod_2(numSeq, residue, n);
                if(!residue) {
                    factors_buf.push_factor(n + 1 - n0, 2);
                }
            }
        }

        constexpr auto max_num_of_primes = size_t{1} << 22;
        constexpr auto segmentLen = uint64_t{1} << 26;

        auto const firstNumber = [&](){
            if constexpr(InitializeFromValue<NumberSequenceT>) {
                return numSeq.value(n0);
            } else {
                return 0;
            }
        }();
        std::vector<uint64_t> primes;
        primes.reserve(max_num_of_primes);

        for (auto fx = f0, fy = (std::min)(f0 + segmentLen, f1);
            fx < f1;
            fx = fy, fy = (std::min)(fy + segmentLen, f1)
        ) {
            primes.clear();
            sieve(fx, fy, primes);

            size_t const iEnd = primes.size();
            #pragma omp parallel for
            for(auto i = 0; i < iEnd; ++i) {
                auto p = primes[i];
                auto barrettMu = compute_barrett_mu_numseq<NumberSequenceT>(p);
                auto residue = [&](){
                    if constexpr(InitializeFromValue<NumberSequenceT>) {
                        return firstNumber.mod(p);
                    } else if constexpr(HasValueModMu<NumberSequenceT>) {
                        return numSeq.value_mod_mu(n0, p, barrettMu);
                    } else {
                        return numSeq.value_mod(n0, p);
                    }
                }();
                if(!residue) {
                    factors_buf.push_factor(0, p);
                }
                /// Propagate residues to next numbers in sequence
                for(auto n = n0, nEnd = n1 - 1; n < nEnd; ++n) {
                    if constexpr(HasNextValueModMu<NumberSequenceT>) {
                        residue = numSeq.next_value_mod_mu(residue, n, p, barrettMu);
                    } else {
                        residue = numSeq.next_value_mod(residue, n, p);
                    }
                    if(!residue) {
                        factors_buf.push_factor(n + 1 - n0, p);
                    }
                }
            }
            updateProgress(fy, !primes.empty() ? primes.back() : 0);
        }
        if(progressHandler) {
            progressHandler->end();
        }
        results = factors_buf.to_factoring_results();
    }


    template <typename NumberSequenceT>
        requires NumberSequence<NumberSequenceT>
    void trial_factor(
        trial_factoring_options const & opts,
        factoring_results<uint64_t, uint32_t> & results,
        std::ostream & out
    )
    {
        if(opts.n1 <= opts.n0) {
            out << "Range [" << opts.n0 << ", " << opts.n1 << "[ is empty, no number to process..." << std::endl;
            return;
        }
        if(opts.f1 <= opts.f0) {
            out << "Range [" << opts.f0 << ", " << opts.f1 << "[ is empty, no prime numbers to divide by..." << std::endl;
            return;
        }

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

        if(results.empty()) {
            results.reserve(opts.n1 - opts.n0);
        } else {
            detail::read_previous_results(opts.n0, opts.n1, results);
        }

        out << "Trial factoring starts, results will be written to "
            << (output_ptr ? (std::string{"file `"} + opts.output_path.value().string() + "'")
                           : "the console.")
            << std::endl;
        out << "TF range is [" << opts.f0 << ", " << opts.f1 << "[." << std::endl;
        out << "Will use up to " << omp_get_max_threads() << " threads on the CPU." << std::endl;

        {
            timer tfTimer{"Trial factoring took ", out};

#ifdef CUTRIALDIVE_ENABLE_GPU
            device_trial_factor<NumberSequenceT>(opts, results, out);
#else
            host_trial_factor<NumberSequenceT>(opts, results, out);
#endif
        }
        (output_ptr ? *output_ptr.get() : out) << results;
    }
}
