/*
* Creation date: 2026.06.03
* Created by Youcef Lemsafer
*/
#pragma once

#include <cstdint>
#include <vector>

namespace cutrialdive {

    enum class PrecomputeReciprocals
    {
        no,
        yes
    };

    template <typename PrimeT, PrecomputeReciprocals havePrecomputedReciprocals>
    class prime_data;

    template <typename PrimeT>
    class prime_data<PrimeT, PrecomputeReciprocals::no>
    {
    public:
        uint64_t size() const;
        void reserve(size_t capacity);
    };
    template <typename PrimeT>
    inline
    uint64_t prime_data<PrimeT, PrecomputeReciprocals::no>::size() const
    {
        return 0;
    }
    template <typename PrimeT>
    inline
    void prime_data<PrimeT, PrecomputeReciprocals::no>::reserve(std::size_t )
    {
    }


    template <typename PrimeT>
    class prime_data<PrimeT, PrecomputeReciprocals::yes>
    {
    public:
        uint64_t size() const;
        void reserve(size_t capacity);
        std::vector<PrimeT> const & reciprocals() const;
        std::vector<PrimeT> & reciprocals();
    private:
        std::vector<PrimeT> reciprocals_;
    };

    template <typename PrimeT>
    inline
    uint64_t prime_data<PrimeT, PrecomputeReciprocals::yes>::size() const
    {
        return reciprocals_.size();
    }

    template <typename PrimeT>
    inline
    void prime_data<PrimeT, PrecomputeReciprocals::yes>::reserve(std::size_t capacity)
    {
        reciprocals_.reserve(capacity);
    }

    template <typename PrimeT>
    inline
    std::vector<PrimeT> const & prime_data<PrimeT, PrecomputeReciprocals::yes>::reciprocals() const
    {
        return reciprocals_;
    }

    template <typename PrimeT>
    inline
    std::vector<PrimeT> & prime_data<PrimeT, PrecomputeReciprocals::yes>::reciprocals()
    {
        return reciprocals_;
    }

    template <typename PrimeT, PrecomputeReciprocals precomputeReciprocals>
    inline void compute_prime_data(std::vector<PrimeT> const & primes, prime_data<PrimeT, precomputeReciprocals> & data)
    {
        if constexpr (precomputeReciprocals == PrecomputeReciprocals::yes) {
            compute_reciprocals(primes, data.reciprocals());
        }
    }

    template <typename PrimeT, PrecomputeReciprocals precomputeReciprocals>
    struct data_impl;

    template <typename PrimeT>
    struct data_impl<PrimeT, PrecomputeReciprocals::yes> {
        using primes_type = PrimeT*;
        using primes_count_type = uint64_t;

        PrimeT * primes;
        PrimeT * reciprocals;
        uint64_t primes_count;
    };

    template <typename PrimeT>
    struct data_impl<PrimeT, PrecomputeReciprocals::no> {
        using primes_type = PrimeT*;
        using primes_count_type = uint64_t;

        primes_type primes;
        primes_count_type primes_count;
    };

    template <typename PrimeT, PrecomputeReciprocals precomputeReciprocals>
    class device_prime_data {
    public:
        device_prime_data() = default;
        device_prime_data(device_prime_data&&) = default;
        device_prime_data(device_prime_data const &) = delete;
        device_prime_data& operator=(device_prime_data const &) = delete;

        data_impl<PrimeT, precomputeReciprocals> const & get_data() const;
        uint64_t size() const;
        void resize(uint64_t newSize);
        void copy_from_host(
                std::vector<PrimeT> const & primes,
                prime_data<PrimeT, precomputeReciprocals> const & data
            );
    private:
        data_impl<PrimeT, precomputeReciprocals> data_{};
    };

    template <typename PrimeT, PrecomputeReciprocals precomputeReciprocals>
    inline
    data_impl<PrimeT, precomputeReciprocals> const & device_prime_data<PrimeT, precomputeReciprocals>::get_data() const
    {
        return data_;   
    }


    template <typename PrimeT, PrecomputeReciprocals precomputeReciprocals>
    inline uint64_t device_prime_data<PrimeT, precomputeReciprocals>::size() const
    {
        return data_.primes_count;
    }

    template <typename PrimeT, PrecomputeReciprocals precomputeReciprocals>
    inline void device_prime_data<PrimeT, precomputeReciprocals>::resize(uint64_t newSize)
    {
        if(newSize <= data_.primes_count) {
            data_.primes_count = newSize;
            return;
        }
        data_.primes_count = newSize;
        free_device_ptr(data_.primes);
        data_.primes = alloc_device_ptr<PrimeT>(data_.primes_count);
        if constexpr (precomputeReciprocals == PrecomputeReciprocals::yes) {
            free_device_ptr(data_.reciprocals);
            data_.reciprocals = alloc_device_ptr<PrimeT>(data_.primes_count);
        }
    }

    template <typename PrimeT, PrecomputeReciprocals precomputeReciprocals>
    inline void device_prime_data<PrimeT, precomputeReciprocals>::copy_from_host(
        std::vector<PrimeT> const & primes,
        prime_data<PrimeT, precomputeReciprocals> const & data
        )
    {
        if constexpr (precomputeReciprocals == PrecomputeReciprocals::yes) {
            // Must have same number of primes and additional data
            // (because there must be a one to one correspondance between
            // i-th prime and associated data)
            assert(primes.size() == data.size());
        }
        auto const primes_size = primes.size();
        resize(primes_size);
        if(!size()) {
            return;
        }
        cudaMemcpy(data_.primes, &primes[0], sizeof(*&primes[0]) * primes_size, cudaMemcpyHostToDevice);
        if constexpr (precomputeReciprocals == PrecomputeReciprocals::yes) {
            cudaMemcpy(data_.reciprocals, &data.reciprocals()[0],
                    sizeof(*&data.reciprocals()[0]) * data.reciprocals().size(),
                    cudaMemcpyHostToDevice
                );
        }
    }
}
