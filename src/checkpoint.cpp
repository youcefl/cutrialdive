/*
* MIT License
* Created on 2026.06.07
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#include "checkpoint.hpp"

#include <utility>
#include <iostream>
#include <fstream>
#include <string>
#include <cassert>

#include "logger.hpp"

namespace {
    using namespace cutrialdive;


    constexpr uint32_t CHECKPOINT_FILE_VERSION = 1;

    template <typename Enum>
    constexpr uint32_t enum_to_u32(Enum e)
    {
        static_assert(std::is_enum_v<Enum>);
        static_assert(
            std::same_as<std::underlying_type_t<Enum>, uint32_t>
        );

        return static_cast<uint32_t>(e);
    }

    template <typename T>
    void write_pod(std::ostream & out, T const & x)
    {
        out.write(reinterpret_cast<char const*>(&x), sizeof(x));
    }

    template <typename T>
    T read_pod(std::istream & in)
    {
        T x;
        in.read(reinterpret_cast<char*>(&x), sizeof(x));
        return x;
    }

    void write_str(std::ostream & out, std::string const & str)
    {
        write_pod(out, str.size());
        out.write(str.data(), str.size());
    }

    std::string read_str(std::istream & in)
    {
        std::string str;
        auto size = read_pod<std::string::size_type>(in);
        str.resize(size);
        in.read(str.data(), size);
        return str;
    }

    void serialize(std::ostream & out, trial_factoring_options const & tfOptions)
    {
        write_pod(out, tfOptions.n0);
        write_pod(out, tfOptions.n1);
        write_pod(out, tfOptions.f0);
        write_pod(out, tfOptions.f1);
        write_str(out, tfOptions.output_path ? tfOptions.output_path->string() : std::string{});
        write_pod(out, tfOptions.max_factors_per_number);
    }

    trial_factoring_options deserialize_tf_options(std::istream & in)
    {
        trial_factoring_options tfOpts;
        tfOpts.n0 = read_pod<decltype(tfOpts.n0)>(in);
        tfOpts.n1 = read_pod<decltype(tfOpts.n1)>(in);
        tfOpts.f0 = read_pod<decltype(tfOpts.f0)>(in);
        tfOpts.f1 = read_pod<decltype(tfOpts.f1)>(in);
        auto outputPath = read_str(in);
        tfOpts.output_path = outputPath.empty() ? std::optional<std::filesystem::path>{}
                                                : std::optional<std::filesystem::path>{outputPath};
        tfOpts.max_factors_per_number = read_pod<decltype(tfOpts.max_factors_per_number)>(in);

        return tfOpts;
    }
}

namespace cutrialdive {

    template <typename PrimeT>
    void serialize_factors_buffer(std::ostream & out, factors_buffer<PrimeT> const & factorsBuf)
    {
        write_pod(out, factorsBuf.numbers_count_);
        write_pod(out, factorsBuf.max_factors_per_number_);
        write_pod(out, factorsBuf.n0_);
        // No need to deal with has_excess_factors_ and last_computed_max_factors_count_
        // they will be reevaluated during TF

        assert(factorsBuf.numbers_count_ == factorsBuf.factors_count_.size());

        for(size_t i = 0, iEnd = factorsBuf.factors_count_.size(); i < iEnd; ++i) {
            auto factorsCount = factorsBuf.factors_count_[i];
            write_pod(out, factorsCount);
            // jEnd is initialized that way because factorsCount is the actual number
            // of factors found which can be greater than max_factors_per_number_
            // (We still count the factors even if we ran out of space to store them)
            for(size_t j = i * factorsBuf.max_factors_per_number_,
                jEnd = j + (std::min)(factorsCount, factorsBuf.max_factors_per_number_);
                j < jEnd;
                ++j) {
                write_pod(out, factorsBuf.factors_[j]);
            }
        }
    }

    template <typename PrimeT>
    factors_buffer<PrimeT> deserialize_factors_buffer(std::istream & in)
    {
        factors_buffer<PrimeT> factorsBuf{0, 0, 0};

        factorsBuf.numbers_count_ = read_pod<decltype(factorsBuf.numbers_count_)>(in);
        factorsBuf.max_factors_per_number_ = read_pod<decltype(factorsBuf.max_factors_per_number_)>(in);
        factorsBuf.n0_ = read_pod<decltype(factorsBuf.n0_)>(in);

        factorsBuf.factors_count_.reserve(factorsBuf.numbers_count_);
        factorsBuf.factors_.resize(factorsBuf.max_factors_per_number_ * factorsBuf.numbers_count_);

        for(size_t i = 0, iEnd = factorsBuf.numbers_count_; i < iEnd; ++i) {
            auto factorsCount = read_pod<typename decltype(factorsBuf.factors_count_)::value_type>(in);
            factorsBuf.factors_count_.push_back(factorsCount);
            // jEnd is initialized that way because factorsCount... etc. (see comment in serialize function above).
            for(size_t j = i * factorsBuf.max_factors_per_number_,
                    jEnd = j + (std::min)(factorsCount, factorsBuf.max_factors_per_number_);
                j < jEnd;
                ++j) {
                auto factor = read_pod<typename decltype(factorsBuf.factors_)::value_type>(in);
                factorsBuf.factors_[j] = factor;
            }
        }

        return std::move(factorsBuf);
    }

    checkpoint_manager::checkpoint_manager(
        std::filesystem::path const & checkpointPath,
        std::chrono::seconds checkpointingPeriod,
        job_spec const & jobSpec
    )
        : checkpoint_path_(checkpointPath)
        , last_write_time_(std::chrono::steady_clock::now())
        , checkpointing_period_(checkpointingPeriod)
        , job_spec_(jobSpec)
    {}

    std::filesystem::path checkpoint_manager::checkpoint_path() const
    {
        return checkpoint_path_;
    }

    std::chrono::seconds checkpoint_manager::period() const
    {
        return checkpointing_period_;
    }

    bool checkpoint_manager::due() const
    {
        return std::chrono::steady_clock::now() - last_write_time_
                        >= checkpointing_period_;
    }

    void checkpoint_manager::write(engine_state engineState)
    {
        auto now = std::chrono::system_clock::now();
        auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(now.time_since_epoch()).count();
        auto tmpChkpntPath = checkpoint_path_;

        tmpChkpntPath += ".tmp_";
        tmpChkpntPath += std::to_string(ns);
        CTDLOG_DEBUG() << "About to write checkpoint to temporary file `" << tmpChkpntPath << "'" << std::endl;
        std::ofstream out{tmpChkpntPath, std::ios::binary};
        if(!out) {
            throw std::runtime_error{"Could not create temporary file `" + tmpChkpntPath.string() + "'"};
        }
        write_pod(out, CHECKPOINT_FILE_VERSION);
        write_pod(out, enum_to_u32(job_spec_.seq_spec.seq_id));
        write_str(out, job_spec_.seq_spec.seq_params);
        serialize(out, job_spec_.tf_options);
        write_pod(out, engineState.last_processed_prime);
        serialize_factors_buffer(out, engineState.factors_buffer.get());
        out.flush();
    
        CTDLOG_DEBUG() << "About to rename file `" << tmpChkpntPath
                       << "' to `" << checkpoint_path_ << "'" << std::endl;
    
        std::filesystem::rename(tmpChkpntPath, checkpoint_path_);
        last_write_time_ = std::chrono::steady_clock::now();
    }

    std::optional<checkpoint_data> checkpoint_manager::load(std::filesystem::path const & checkpointPath)
    {
        CTDLOG_DEBUG() << "About to read checkpoint file at `" << checkpointPath << "'" << std::endl;

        std::ifstream in{checkpointPath, std::ios::binary};
        if(!in) {
            CTDLOG_DEBUG() << "Could not read checkpoint file at `" << checkpointPath << "'" << std::endl;
            return {};
        }
        auto version = read_pod<uint32_t>(in);
        if(version != CHECKPOINT_FILE_VERSION) {
            throw std::runtime_error{"Unexpected version number `" + std::to_string(version) + "' found in checkpoint file"};
        }
        auto numSeqSpec = num_seq_spec {
            static_cast<num_seq_id>(read_pod<uint32_t>(in)),
            read_str(in)
        };
        auto tfOpts = deserialize_tf_options(in);
        auto lastProcessedPrime = read_pod<decltype(std::declval<engine_state>().last_processed_prime)>(in);
        auto factorsBuf = deserialize_factors_buffer<uint64_t>(in);

        return checkpoint_data{
            job_spec{numSeqSpec, tfOpts},
            lastProcessedPrime, factors_buffer_holder{std::move(factorsBuf)}
        };
    }

    void checkpoint_manager::remove_checkpoint()
    {
        CTDLOG_DEBUG() << "Removing file `" << checkpoint_path_ << "'" << std::endl;
        std::filesystem::remove(checkpoint_path_);
    }
}

