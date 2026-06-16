/*
* Youcef Lemsafer, 2026.04.21
* Definition of the sieving back-end
*/
#pragma once

#include <cstdint>
#include <vector>
#include <span>
#ifdef CUTRIALDIVE_HAS_PRIMESIEVE
#include <vector>
#include <thread>
#include <mutex>
#include <memory>
#include <algorithm>
#include <condition_variable>
#include <primesieve.hpp>
#endif
#ifdef CUTRIALDIVE_HAS_LFP
#include <lfp.hpp>
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

    template <sieve_lib SieveLib> class sieve_backend;

    template <sieve_lib SieveLib>
    class siever_impl
    {
    public:
        siever_impl(uint32_t threadsCount);

        /// Sieves range [start, end[ for primes
        /// @param[in] start range start
        /// @param[in] end range end
        /// @param[out] primes container receiving the primes
        void sieve(uint64_t start, uint64_t end, std::vector<std::span<uint64_t>> & primes);

        uint32_t threads_count() const;

    private:
        sieve_backend<SieveLib> siever_;
    };

    using siever = siever_impl<
#ifdef CUTRIALDIVE_HAS_PRIMESIEVE
                        primesieve
#elif defined CUTRIALDIVE_HAS_LFP
                        lfp
#endif
        >;
}


namespace cutrialdive {

    namespace details {

        inline
        uint32_t safe_threads_count(uint32_t threadsCount)
        {
            auto result = threadsCount ? threadsCount : std::thread::hardware_concurrency();
            return result ? result : 1;
        }

    }

#ifdef CUTRIALDIVE_HAS_PRIMESIEVE

    template <>
    class sieve_backend<sieve_lib::primesieve>
    {
    public:
        sieve_backend(sieve_backend<sieve_lib::primesieve>&&) noexcept = default;
        sieve_backend<sieve_lib::primesieve>& operator=(sieve_backend<sieve_lib::primesieve>&&) noexcept = default;
        sieve_backend(sieve_backend<sieve_lib::primesieve> const &) = delete;
        sieve_backend<sieve_lib::primesieve>& operator=(sieve_backend<sieve_lib::primesieve> const &) = delete;
        explicit sieve_backend(uint32_t threadsCount);
        ~sieve_backend();
        void sieve(uint64_t start, uint64_t end, std::vector<std::span<uint64_t>> & primes);
        uint32_t threads_count() const;

    private:
        struct range
        {
            uint64_t start;
            uint64_t end;
        };
        struct sieved_range
        {
            range input_range;
            std::vector<uint64_t> primes;
        };
        uint32_t threads_count_;
        std::vector<std::thread> sievers_;
        std::vector<range> ranges_;
        std::mutex ranges_mutex_;
        std::vector<sieved_range> sieved_ranges_;
        std::atomic_uint32_t sieved_ranges_count_;
        std::condition_variable work_added_;
        std::condition_variable work_done_;
        uint32_t pending_ranges_count_;
        std::mutex pending_ranges_mutex_;
        bool is_terminating_;
    };

    inline
    sieve_backend<sieve_lib::primesieve>::sieve_backend(uint32_t threadsCount)
        : threads_count_(details::safe_threads_count(threadsCount))
        , sieved_ranges_count_(0)
        , pending_ranges_count_(0)
        , is_terminating_(false)
    {
        sieved_ranges_.resize(threads_count_);
        std::for_each(std::begin(sieved_ranges_), std::end(sieved_ranges_), [this](auto & sievedRange) {
            sievedRange.primes.reserve((1u << 22) / threads_count_);
        });
        if(threads_count_ == 1) {
            // In this case use main thread, do not create new threads
            return;
        }
        for(uint32_t i = 0; i < threads_count_; ++i) {
            sievers_.emplace_back(std::thread{[this](){
                for(;;) {
                    std::unique_lock lck{ranges_mutex_};
                    work_added_.wait(lck, [&]{ return !ranges_.empty() || is_terminating_; });
                    if(is_terminating_) {
                        return;
                    }
                    auto range = ranges_.back();
                    ranges_.pop_back();
                    lck.unlock();
                    {
                        auto & sievedRange = sieved_ranges_[
                                                sieved_ranges_count_.fetch_add(1, std::memory_order_relaxed)];
                        sievedRange.input_range.start = range.start;
                        sievedRange.input_range.end = range.end;
                        sievedRange.primes.clear();
                        primesieve::generate_primes(range.start, range.end, &sievedRange.primes);
                    }
                    {
                        std::unique_lock lck{pending_ranges_mutex_};
                        --pending_ranges_count_;
                    }
                    work_done_.notify_one();
                }
            }});
        }
    }

    inline
    uint32_t sieve_backend<sieve_lib::primesieve>::threads_count() const
    {
        return threads_count_;
    }

    inline
    sieve_backend<sieve_lib::primesieve>::~sieve_backend()
    {
        is_terminating_ = true;
        work_added_.notify_all();
        std::for_each(std::begin(sievers_), std::end(sievers_), [](auto & siever) {
            siever.join();
        });
    }

    inline
    void sieve_backend<sieve_lib::primesieve>::sieve(
        uint64_t start,
        uint64_t end,
        std::vector<std::span<uint64_t>> & primes
      )
    {
        primes.clear();
        if(start >= end) {
            return;
        }
        if(threads_count_ == 1) {
            sieved_ranges_[0].input_range.start = start;
            sieved_ranges_[0].input_range.end = end;
            sieved_ranges_[0].primes.clear();
            primesieve::generate_primes(start, end, &sieved_ranges_[0].primes);
            if(!sieved_ranges_[0].primes.empty()) {
                primes.emplace_back(&sieved_ranges_[0].primes[0], sieved_ranges_[0].primes.size());
            }
            return;
        }
        auto const rangeBaseLength = (end - start) >= 2 * threads_count_ ? (end - start) / threads_count_ : end - start;
        auto const rangeCorrection = (end - start) >= 2 * threads_count_ ? (end - start) % threads_count_ : 0;
        auto const chunksCount = (end - start) >= 2 * threads_count_ ? threads_count_ : 1;

        // This lock is here to make it impossible for the workers to notify end of work
        // before main thread is ready to wait (i.e. all ranges are pushed).
        std::unique_lock lckPendingRanges{pending_ranges_mutex_};
        pending_ranges_count_ = chunksCount;
        sieved_ranges_count_.store(0, std::memory_order_relaxed);
        if(sieved_ranges_.size() != chunksCount) {
            sieved_ranges_.resize(chunksCount);
        }
        auto rangeIndex  = 0;
        for(auto rangeStart = start, rangeLength = rangeBaseLength + (rangeCorrection ? 1 : 0);
           rangeStart < end;
           rangeStart = (std::min)(rangeStart + rangeLength, end),
           rangeLength = (rangeIndex < rangeCorrection) ? rangeBaseLength + 1 : rangeBaseLength,
           rangeIndex++
        ) {
            {
                std::lock_guard<std::mutex> lck(ranges_mutex_);
                ranges_.push_back(range{rangeStart, (rangeIndex < chunksCount - 1) ? rangeStart + rangeLength - 1 : end - 1});
            }
            work_added_.notify_one();
        }

        work_done_.wait(lckPendingRanges, [this]{ return !pending_ranges_count_; });
        std::sort(std::begin(sieved_ranges_), std::begin(sieved_ranges_) + sieved_ranges_count_,
            [](auto & rx, auto & ry) {
                return rx.input_range.start < ry.input_range.start;
            });
        for(auto i = 0u; i < sieved_ranges_count_; ++i) {
            if(sieved_ranges_[i].primes.empty()) {
                continue;
            }
            primes.emplace_back(&sieved_ranges_[i].primes[0], sieved_ranges_[i].primes.size());
        }
    }

#endif // CUTRIALDIVE_HAS_PRIMESIEVE

#ifdef CUTRIALDIVE_HAS_LFP

    template <>
    class sieve_backend<sieve_lib::lfp>
    {
    public:
        sieve_backend(uint32_t threadsCount);
        void sieve(uint64_t start, uint64_t end, std::vector<std::span<uint64_t>> & primes);
        uint32_t threads_count() const;
    private:
        uint32_t threads_count_;
        std::vector<uint64_t> primes_;
    };

    inline
    sieve_backend<sieve_lib::lfp>::sieve_backend(uint32_t threadsCount)
        : threads_count_(details::safe_threads_count(threadsCount))
    {}

    inline
    uint32_t sieve_backend<sieve_lib::lfp>::threads_count() const
    {
        return threads_count_;
    }

    inline
    void sieve_backend<sieve_lib::lfp>::sieve(uint64_t start, uint64_t end, std::vector<std::span<uint64_t>> & primes)
    {
        primes.clear():
        auto sieveRes = lfp::sieve<uint64_t>(start, end, lfp::threads{threads_count_});
        auto rng = sieveRes.range();
        primes_.assign(rng.begin(), rng.end());
        primes.emplace_back(&primes_[0], primes_.size());
    }

#endif // CUTRIALDIVE_HAS_LFP

    template <sieve_lib SieveLib>
    inline
    siever_impl<SieveLib>::siever_impl(uint32_t threadsCount)
        : siever_(threadsCount)
    {}

    template <sieve_lib SieveLib>
    inline
    void siever_impl<SieveLib>::sieve(uint64_t start, uint64_t end, std::vector<std::span<uint64_t>> & primes)
    {
        siever_.sieve(start, end, primes);
    }

    template <sieve_lib SieveLib>
    inline
    uint32_t siever_impl<SieveLib>::threads_count() const
    {
        return siever_.threads_count();
    }
}
