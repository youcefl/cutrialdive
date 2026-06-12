/*
 * Youcef Lemsafer
 * 2023.01.19
 */
#include "prp.hpp"

#include <iostream>

#include "hgint.hpp"
#include "gw_utility.h"
#include "number_sequence.hpp"
#include "builtin_number_sequences.hpp"
#include "number_sequence_helpers.hpp"
#include "timer.hpp"
#include "num_seq_dispatch.hpp"

bool isPRP(HgInt number)
{
    return gw_prp(number.get()) != 0;
}

bool is_prp(HgInt number)
{
    return isPRP(number);
}

namespace cutrialdive {

    template <typename NumberSequenceT>
        requires NumberSequence<NumberSequenceT>
    void run_prp_test(uint64_t n,
        std::vector<factor<uint64_t, uint32_t>> & factors,
        bool haveToBoostFactors,
        std::ostream & out,
        auto&&... args
    )
    {
        NumberSequenceT numSeq{std::forward<decltype(args)>(args)...};
        std::string const numberAsStr = short_name(numSeq) + ("(" + std::to_string(n)) + ")";
        auto pureSeq = get_math_sequence(numSeq);
        HgInt cofactor;
        {
            auto number = pureSeq.value(n);
            timer t{(std::string{"Computing the cofactor of "}
                        + numberAsStr + " took ").c_str(), out};
            out << "Boosting factors: " << std::boolalpha << haveToBoostFactors << std::endl;
            cofactor = haveToBoostFactors ? compute_cofactor_boosted(number, factors)
                        : compute_cofactor_exact(number, factors);
            out << "Factors: ";
            char const * sep = "";
            std::for_each(std::begin(factors), std::end(factors), [&out, &sep](auto const & fact) {
                out << sep << fact;
                sep = ", ";
            });
            out << std::endl;
        }
        {
            auto sizeInBase10 = cofactor.sizeInBase(10);
            timer t{(std::string{"The PRP test on the " + std::to_string(sizeInBase10) 
                + " digits cofactor of "} + numberAsStr
                + " took ").c_str(), out};
            auto isPrp = is_prp(std::move(cofactor));
            out << "Is PRP: " << (isPrp ? "YES" : "NO") << std::endl;
        }
    }


    void run_prp_test(
        num_seq_spec numSeqSpec,
        uint64_t n,
        std::vector<factor<uint64_t, uint32_t>> & factors,
        bool haveToBoostFactors,
        std::ostream & out
    )
    {
        timer tfTimer{"PRP test took ", std::cout};

        dispatch_num_seq<decltype(n)>(numSeqSpec, [&]<typename Seq, typename... Args>(Args&&... args) {
            run_prp_test<Seq>(n, factors, haveToBoostFactors, std::cout, std::forward<Args>(args)...);
        });
    }

}
