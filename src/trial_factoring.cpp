/*
* Created on 2026.05.04
* Copyright (c) Youcef Lemsafer
*/
#include "trial_factoring.hpp"
#include "trial_factoring_impl.hpp"
#include "checkpoint.hpp"


namespace fs = std::filesystem;

namespace cutrialdive {

    namespace {

        int const defaultCheckpointingPeriodSeconds = 360;

        factoring_results<uint64_t, uint32_t>
        trial_factor(
            num_seq_spec numSeqSpec,
            trial_factoring_options const & opts,
            std::ostream & out,
            checkpoint_manager * checkpoint,
            engine_state * resumeState = nullptr
        )
        {
            factoring_results<uint64_t, uint32_t> results{opts.n0, opts.n1 - opts.n0};
            dispatch_num_seq<decltype(opts.n0)>(numSeqSpec, [&]<typename Seq>() {
                trial_factor<Seq>(opts, results, out, checkpoint, resumeState);
            });
            return results;
        }

    }

    /// @brief Performs trial factoring according to given options
    /// @param opts trial factoring options
    factoring_results<uint64_t, uint32_t>
    trial_factor(
        num_seq_spec numSeqSpec,
        trial_factoring_options const & opts,
        std::ostream & out
    )
    {
        std::unique_ptr<checkpoint_manager> checkpoint;
        if(opts.output_path) {
            checkpoint = std::make_unique<checkpoint_manager>(
                std::filesystem::path{opts.output_path.value()} += ".chkpnt",
                std::chrono::seconds{defaultCheckpointingPeriodSeconds},
                job_spec{numSeqSpec, opts}
            );
        }
        return trial_factor(numSeqSpec, opts, out, checkpoint.get());
    }

    /// Resumes an interrupted trial factoring
    /// @param[in] checkpointPath path of the checkpoint file
    /// @param[in/out] out stream where to write user messages
    factoring_results<uint64_t, uint32_t>
    resume_trial_factoring(
        std::filesystem::path const & checkpointPath,
        std::ostream & out
    )
    {
        out << "Resuming from checkpoint `" << checkpointPath.string() << "'..." << std::endl;
        auto checkpointData = checkpoint_manager::load(checkpointPath);
        if(!checkpointData) {
            throw std::runtime_error{std::string{"Failed to load checkpoint data from file `"}
                                        +  checkpointPath.string() + "'"};
        }
        auto chkpntPath = fs::path{checkpointData->spec.tf_options.output_path.value()} += ".chkpnt";
        checkpoint_manager checkpoint{chkpntPath,
                std::chrono::seconds{defaultCheckpointingPeriodSeconds},
                checkpointData->spec
            };

        return trial_factor(
                checkpointData->spec.seq_spec,
                checkpointData->spec.tf_options,
                out,
                &checkpoint,
                &checkpointData->state
            );
    }
}
