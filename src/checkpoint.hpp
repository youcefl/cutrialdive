/*
* Created on 2026.06.07
* Copyright (c) Youcef Lemsafer
*/
#pragma once

#include <cstdint>
#include <chrono>
#include <optional>
#include <filesystem>
#include <memory>

#include "num_seq_spec.hpp"
#include "trial_factoring_options.hpp"
#include "factors_buffer.hpp"


namespace cutrialdive {

    /// Holds a factors buffer regardless of ownership
    class factors_buffer_holder
    {
    public:
        factors_buffer_holder(factors_buffer_holder const &) = delete;
        factors_buffer_holder& operator=(factors_buffer_holder const &) = delete;
        factors_buffer_holder(factors_buffer_holder&&) noexcept;
        factors_buffer_holder& operator=(factors_buffer_holder&&) noexcept;

        /// Constructs an instance with a non-owned factors buffer
        explicit factors_buffer_holder(factors_buffer<uint64_t> * factorsBuffer);
        /// Constructs an instance with an owned factors buffer
        explicit factors_buffer_holder(factors_buffer<uint64_t> factorsBuffer);
        /// Returns the underlying factors buffer as a const
        factors_buffer<uint64_t> const & get() const;
        /// Returns the underlying factors buffer
        factors_buffer<uint64_t> & get();
    private:
        factors_buffer<uint64_t> * factors_buffer_;
        std::optional<factors_buffer<uint64_t>> owned_factors_buffer_;
    };

    struct job_spec
    {
        /// Number sequence specification
        num_seq_spec seq_spec;
        /// Trial factoring options
        trial_factoring_options tf_options;
    };

    struct engine_state
    {
        /// The last prime successfully processed
        uint64_t last_processed_prime;
        /// Factors found during processing of primes <= last_processed_prime.
        factors_buffer_holder factors_buffer;
    };

    struct checkpoint_data
    {
        job_spec spec;
        engine_state state;
    };

    class checkpoint_manager
    {
    public:
        /// Constructs an instance with a path, a checkpointing period and a job specification
        /// @param[in] checkpointPath path of the checkpoint file
        /// @param[in] checkpointingPeriod checkpointing period in seconds
        /// @param[in] jobSpec job specification
        checkpoint_manager(
            std::filesystem::path const & checkpointPath,
            std::chrono::seconds checkpointingPeriod,
            job_spec const & jobSpec
        );

        /// Returns the path of the checkpoint file
        std::filesystem::path checkpoint_path() const;

        /// Returns true if (and only if) a checkpoint file has to be written
        bool due() const;

        /// Writes a checkpoint @c checkpoint_path()
        void write(engine_state engineState);

        /// Loads a checkpoint from @param checkpointPath
        static std::optional<checkpoint_data> load(std::filesystem::path const & checkpointPath);

        void end();
    private:
        job_spec job_spec_;
        std::filesystem::path checkpoint_path_;
        std::chrono::time_point<std::chrono::steady_clock> last_write_time_;
        std::chrono::seconds checkpointing_period_;
    };
}


namespace cutrialdive {

    inline
    factors_buffer_holder::factors_buffer_holder(factors_buffer<uint64_t> * factorsBuffer)
        : factors_buffer_(factorsBuffer)
    {}

    inline
    factors_buffer_holder::factors_buffer_holder(factors_buffer<uint64_t> factorsBuffer)
        : owned_factors_buffer_(std::move(factorsBuffer))
    {
        factors_buffer_ = &owned_factors_buffer_.value();
    }

    inline
    factors_buffer_holder::factors_buffer_holder(factors_buffer_holder&& other) noexcept
        : owned_factors_buffer_(std::move(other.owned_factors_buffer_))
    {
        factors_buffer_ = owned_factors_buffer_ ? &*owned_factors_buffer_ : other.factors_buffer_;
    }

    inline
    factors_buffer_holder&
    factors_buffer_holder::operator =(factors_buffer_holder&& other) noexcept
    {
        owned_factors_buffer_ = std::move(other.owned_factors_buffer_);
        factors_buffer_ = owned_factors_buffer_ ? &*owned_factors_buffer_ : other.factors_buffer_;
        return *this;
    }

    inline
    factors_buffer<uint64_t> const & factors_buffer_holder::get() const
    {
        return *factors_buffer_;
    }

    inline
    factors_buffer<uint64_t> & factors_buffer_holder::get()
    {
        return *factors_buffer_;
    }

}