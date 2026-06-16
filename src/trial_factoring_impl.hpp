/*
* Created on 2026.06.04
* Copyright (c) Youcef Lemsafer
*/
#pragma once

#include <cstdint>
#include <type_traits>
#include <set>
#include <memory>
#include <atomic>
#include <mutex>
#include <condition_variable>
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
#include "checkpoint.hpp"
#include "trial_factoring_context.hpp"


namespace cutrialdive {
namespace details {

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

    struct chunk_info
    {
        uint64_t first_prime;
        uint64_t previous_chunk_last_prime;
    };

    inline bool operator <(chunk_info const & x, chunk_info const & y)
    {
        return x.first_prime < y.first_prime;
    }

} // namespace cutrialdive::detail


    template <typename NumberSequenceT>
    class trial_factorer
    {
    public:
        trial_factorer(trial_factorer const &) = delete;
        trial_factorer& operator=(trial_factorer const &) = delete;

        trial_factorer(
            trial_factoring_context& ctx,
            NumberSequenceT & numSeq,
            factors_buffer<uint64_t> & factorsBuf,
            std::function<void(uint64_t, uint64_t)> updateProgress
        );
        ~trial_factorer();
        void process_chunks(std::vector<std::span<uint64_t>> primes);

    private:
        void process_chunk(std::span<uint64_t> primes);
        static auto first_value(NumberSequenceT & numSeq, typename NumberSequenceT::index_type n);
        

        uint32_t threads_count_;
        typename NumberSequenceT::index_type n0_;
        NumberSequenceT & num_seq_;
        using first_value_type = decltype(first_value(num_seq_, typename NumberSequenceT::index_type{}));
        first_value_type first_number_;
        factors_buffer<uint64_t> & factors_buf_;
        trial_factoring_context& context_;
        std::vector<std::thread> chunk_processors_;
        std::deque<std::span<uint64_t>> chunk_queue_;
        std::condition_variable chunk_enqueued_;
        std::condition_variable chunk_processed_;
        std::mutex chunk_queue_mutex_;
        std::set<details::chunk_info> pending_chunks_;
        std::mutex progress_mutex_;
        std::atomic_size_t pending_chunks_count_;
        std::atomic_bool no_more_work_;
        std::function<void(uint64_t, uint64_t)> update_progress_;
    };

    template <typename NumberSequenceT>
    inline
    trial_factorer<NumberSequenceT>::trial_factorer(
        trial_factoring_context& ctx,
        NumberSequenceT & numSeq,
        factors_buffer<uint64_t> & factorsBuf,
        std::function<void(uint64_t, uint64_t)> updateProgress
      ) : threads_count_(ctx.runtime_options.threads_count)
        , n0_(ctx.options.n0)
        , num_seq_(numSeq)
        , first_number_(first_value(num_seq_, n0_))
        , factors_buf_(factorsBuf)
        , context_(ctx)
        , pending_chunks_count_{}
        , no_more_work_{}
        , update_progress_(updateProgress)
    {
        if(threads_count_ == 1) {
            // Use main thread
            return;
        }
        for(auto i = 0u; i < threads_count_; ++i) {
            chunk_processors_.emplace_back(std::thread{[this](){
                for(;;) {
                    std::unique_lock<std::mutex> lck{chunk_queue_mutex_};
                    chunk_enqueued_.wait(lck, [&]{ return no_more_work_.load(std::memory_order_relaxed)
                                                        || !chunk_queue_.empty(); });
                    if(no_more_work_.load(std::memory_order_relaxed)) {
                        return;
                    }
                    auto chunk = chunk_queue_.front();
                    chunk_queue_.pop_front();
                    lck.unlock();
                    process_chunk(chunk);
                    {
                        std::unique_lock<std::mutex> lckPendingChunks{progress_mutex_};
                        pending_chunks_.erase({chunk.front(), 0});
                        if(!pending_chunks_.empty()) {
                            auto firstPendingChunkInfo = *std::begin(pending_chunks_);
                            if(firstPendingChunkInfo.previous_chunk_last_prime) {
                                update_progress_(firstPendingChunkInfo.first_prime,
                                    firstPendingChunkInfo.previous_chunk_last_prime);
                            }
                        }
                    }
                    pending_chunks_count_.fetch_sub(1, std::memory_order_relaxed);
                    chunk_processed_.notify_one();
                }
            }});
        }
    }

    template <typename NumberSequenceT>
    inline
    trial_factorer<NumberSequenceT>::~trial_factorer()
    {
        if(!chunk_processors_.empty()) {
            no_more_work_.store(true, std::memory_order_relaxed);
            chunk_enqueued_.notify_all();
            std::for_each(std::begin(chunk_processors_), std::end(chunk_processors_), [](auto & t) { t.join(); });
        }
    }

    template <typename NumberSequenceT>
    inline
    auto trial_factorer<NumberSequenceT>::first_value(
        NumberSequenceT & numSeq,
        typename NumberSequenceT::index_type n
      )
    {
        if constexpr(InitializeFromValue<NumberSequenceT>) {
            return numSeq.value(n);
        } else {
            return 0;
        }
    }

    template <typename NumberSequenceT>
    inline
    void trial_factorer<NumberSequenceT>::process_chunk(std::span<uint64_t> primes)
    {
        auto const n1 = context_.options.n1;
        size_t const iEnd = primes.size();
        for(auto i = 0; i < iEnd; ++i) {
            auto p = primes[i];
            auto barrettMu = compute_barrett_mu_numseq<NumberSequenceT>(p);
            auto residue = [&](){
                if constexpr(InitializeFromValue<NumberSequenceT>) {
                    return first_number_.mod(p);
                } else if constexpr(HasValueModMu<NumberSequenceT>) {
                    return num_seq_.value_mod_mu(n0_, p, barrettMu);
                } else {
                    return num_seq_.value_mod(n0_, p);
                }
            }();
            if(!residue) {
                factors_buf_.push_factor(0, p);
            }
            /// Propagate residues to next numbers in sequence
            for(auto n = n0_, nEnd = n1 - 1; n < nEnd; ++n) {
                if constexpr(HasNextValueModMu<NumberSequenceT>) {
                    residue = num_seq_.next_value_mod_mu(residue, n, p, barrettMu);
                } else {
                    residue = num_seq_.next_value_mod(residue, n, p);
                }
                if(!residue) {
                    factors_buf_.push_factor(n + 1 - n0_, p);
                }
            }
        }
    }

    template <typename NumberSequenceT>
    inline
    void trial_factorer<NumberSequenceT>::process_chunks(std::vector<std::span<uint64_t>> primes)
    {
        if(threads_count_ == 1) {
            for(auto chunk : primes) {
                process_chunk(chunk);
            }
            return;
        }
        pending_chunks_count_.store(primes.size(), std::memory_order_relaxed);
        pending_chunks_.clear();
        auto prevChunkLastPrime = uint64_t{};
        for(auto chunk : primes) {
            pending_chunks_.insert({chunk.front(), prevChunkLastPrime});
            prevChunkLastPrime = chunk.back();
        }
        for(auto chunk : primes) {
            std::unique_lock<std::mutex> lck{chunk_queue_mutex_};
            chunk_queue_.push_back(chunk);
            lck.unlock();
            chunk_enqueued_.notify_one();
        }
        // Wait until all workers are done
        std::unique_lock<std::mutex> lck{chunk_queue_mutex_};
        chunk_processed_.wait(lck, [&]{ return !pending_chunks_count_.load(std::memory_order_relaxed); });
    }

    template <typename NumberSequenceT>
        requires PureMathSequence<NumberSequenceT>
    void host_trial_factor(
        trial_factoring_context& ctx,
        NumberSequenceT numSeq
        )
    {
        auto const & opts = ctx.options;
        auto & out = ctx.output_stream;
        auto resumeState = ctx.resume_state;
        auto & results = ctx.results;
        uint64_t n0 = opts.n0, n1 = opts.n1;
        uint64_t f0 = resumeState ? resumeState->last_processed_prime + 1 : opts.f0;
        uint64_t f1 = opts.f1;

        if(resumeState) {
            out << "Resuming at the smallest prime > " << resumeState->last_processed_prime << "." << std::endl;
        }

        factors_buffer<uint64_t> factorsBuf = resumeState 
                ? std::move(resumeState->factors_buffer.get())
                : factors_buffer<uint64_t>{n0, n1 - n0, opts.max_factors_per_number};

        auto progressHandler = opts.is_progress_enabled 
                ? std::make_unique<progress>(f1, ctx.runtime_options.progress_period, out)
                : std::unique_ptr<progress>{};
        auto updateProgress = progressHandler ? std::function<void(uint64_t, uint64_t)>{
                    [&progressHandler](auto segUpperBound, auto lastPrimeInSeg) {
                        progressHandler->update(segUpperBound, lastPrimeInSeg); 
                    }}
                : std::function<void(uint64_t, uint64_t)>{[](auto, auto) {}};

        // 2 is a special case, deal with it now so that only odd primes remain
        if(f0 <= 2 && f1 > 2) {
            f0 = 3;
            auto residue = value_mod_2(numSeq, n0);
            if(!residue) {
                factorsBuf.push_factor(0, 2);
            }
            for(auto n = n0, nEnd = n1 - 1; n < nEnd; ++n) {
                residue = next_value_mod_2(numSeq, residue, n);
                if(!residue) {
                    factorsBuf.push_factor(n + 1 - n0, 2);
                }
            }
        }

        constexpr auto segmentLen = uint64_t{1} << 26;
        std::vector<std::span<uint64_t>> primes;
        auto checkpoint = ctx.checkpoint;
        siever primeGen{ctx.runtime_options.threads_count};

        trial_factorer<NumberSequenceT> factorer{ctx, numSeq, factorsBuf, updateProgress};

        for (auto fx = f0, fy = (std::min)(f0 + segmentLen, f1);
            fx < f1;
            fx = fy, fy = (std::min)(fy + segmentLen, f1)
        ) {
            primes.clear();
            primeGen.sieve(fx, fy, primes);

            factorer.process_chunks(primes);

            auto lastPrime = !primes.empty() ? primes.back().back() : 0;
            updateProgress(fy, lastPrime);
            if(!primes.empty() && checkpoint && checkpoint->due()) {
                checkpoint->write(engine_state{
                    lastPrime,
                    factors_buffer_holder{&factorsBuf}
                });
            }
        }
        if(progressHandler) {
            progressHandler->end();
        }
        if(checkpoint) {
            checkpoint->end();
        }
        ctx.results = factorsBuf.to_factoring_results<uint64_t, uint32_t>();
    }


    template <typename NumberSequenceT>
        requires NumberSequence<NumberSequenceT>
    void trial_factor(
        trial_factoring_context&& ctx,
        auto&&... args
    )
    {
        auto const & opts = ctx.options;
        auto & out = ctx.output_stream;
        auto & results = ctx.results;

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
            if(!ctx.checkpoint && std::filesystem::exists(path)) {
                std::string msg{"Output: cannot overwrite existing file "};
                msg += path.string();
                throw std::runtime_error(msg);
            }
            output_ptr = std::make_unique<std::ofstream>(path);
        }

        if(results.empty()) {
            results.reserve(opts.n1 - opts.n0);
        } else {
            details::read_previous_results(opts.n0, opts.n1, results);
        }

        out << "Trial factoring starts, results will be written to "
            << (output_ptr ? (std::string{"file `"} + opts.output_path.value().string() + "'")
                           : "the console.")
            << std::endl;
        out << "TF range is [" << opts.f0 << ", " << opts.f1 << "[." << std::endl;
        out << "Will use up to " << omp_get_max_threads() << " threads on the CPU." << std::endl;
        auto checkpoint = ctx.checkpoint;
        if(checkpoint) {
            auto checkpointPath = checkpoint->checkpoint_path();
            out << "Checkpoint file will be written to `" << checkpointPath.string() 
                << "' (period is " << checkpoint->period() << ")." << std::endl;
            out << "Use `--resume " << checkpointPath << "' to resume execution." << std::endl;
        }
        {
            timer tfTimer{"Trial factoring took ", out};

            NumberSequenceT numSeq{std::forward<decltype(args)>(args)...};
            auto pureSeq = get_math_sequence(numSeq);
#ifdef CUTRIALDIVE_ENABLE_GPU
            device_trial_factor(ctx, pureSeq);
#else
            host_trial_factor(ctx, pureSeq);
#endif
        }
        (output_ptr ? *output_ptr.get() : out) << results;
    }
}
