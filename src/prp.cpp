/*
 * Youcef Lemsafer
 * 2023.01.19
 */
#include "prp.hpp"
#include "hgint.hpp"
#include "gw_utility.h"
#include "number_sequence.hpp"
#include "builtin_number_sequences.hpp"
#include "timer.hpp"
#include "mode.hpp"

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
        std::ostream & out
    )
    {
        std::string const numberAsStr = NumberSequenceT::short_name() + ("(" + std::to_string(n)) + ")";
        HgInt cofactor;
        {
            auto number = NumberSequenceT::value(n);
            timer t{(std::string{"Computing the cofactor of "}
                        + numberAsStr + " took ").c_str(), out};
            out << "Boosting factors: " << std::boolalpha << haveToBoostFactors << std::endl;
            out << "Factors: ";
            std::for_each(std::begin(factors), std::end(factors), [&out](auto const & fact) {
                out << ", " << fact;
            });
            out << std::endl;
            cofactor = haveToBoostFactors ? compute_cofactor_boosted(number, factors)
                        : compute_cofactor_exact(number, factors);
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
        mode_flag modeFlag,
        uint64_t n,
        std::vector<factor<uint64_t, uint32_t>> & factors,
        bool haveToBoostFactors,
        std::ostream & out
    )
    {
        timer tfTimer{"PRP test took ", std::cout};

        dispatch_mode<decltype(n)>(modeFlag, [&]<typename Seq>() {
            run_prp_test<Seq>(n, factors, haveToBoostFactors, std::cout);
        });
    }

}
