/*
* Created on 2026.06.09
* Copyright (c) Youcef Lemsafer
*/
#include <filesystem>
#include <iostream>
#include <chrono>

#include <catch2/catch_test_macros.hpp>

#include "checkpoint.hpp"
#include "num_seq_spec.hpp"
#include "trial_factoring_options.hpp"
#include "factors_buffer.hpp"
#include "filesystem_helper.hpp"

using namespace cutrialdive;
namespace fs = std::filesystem;


TEST_CASE("Basic checkpoint test")
{
    constexpr auto maxFactorsPerNumber = 32;

    num_seq_spec numSeqSpec{num_seq_id::mersenne};
    trial_factoring_options tfOpts{101, 105, 0, 1u << 16, {},
            maxFactorsPerNumber};

    auto checkpointPath = tests::get_temp_file_path() += ".chkpnt";

    checkpoint_manager checkpoint{checkpointPath, std::chrono::seconds{1800}, job_spec{numSeqSpec, tfOpts}};
    
    factors_buffer<uint64_t> factorsBuf{tfOpts.n0, tfOpts.n1 - tfOpts.n0, tfOpts.max_factors_per_number};

    factorsBuf.push_factor(102 - 101, 3);
    factorsBuf.push_factor(102 - 101, 7);
    factorsBuf.push_factor(102 - 101, 103);
    factorsBuf.push_factor(104 - 101, 3);
    factorsBuf.push_factor(104 - 101, 5);
    factorsBuf.push_factor(104 - 101, 17);
    factorsBuf.push_factor(104 - 101, 53);

    constexpr auto lastProcessedPrime = 127;

    checkpoint.write(engine_state{
        .last_processed_prime = lastProcessedPrime,
        .factors_buffer = factors_buffer_holder{&factorsBuf}
    });

    auto optLoaded = checkpoint_manager::load(checkpointPath);

    REQUIRE(optLoaded.has_value());

    auto const & loaded = *optLoaded;

    REQUIRE(loaded.state.last_processed_prime == lastProcessedPrime);
    REQUIRE(loaded.spec.seq_spec == numSeqSpec);
    REQUIRE(loaded.spec.tf_options == tfOpts);
    REQUIRE(loaded.state.factors_buffer.get().n0() == 101);
    REQUIRE(loaded.state.factors_buffer.get().numbers_count() == 4);
    REQUIRE(loaded.state.factors_buffer.get().max_factors_per_number() == maxFactorsPerNumber);

    auto factoringResults = loaded.state.factors_buffer.get().to_factoring_results<uint64_t, uint32_t>();

    REQUIRE(factoringResults.size() == 4);
    REQUIRE(factoringResults[0].size() == 0);

    checkpoint.remove_checkpoint();

    REQUIRE(!fs::exists(checkpointPath));
}
