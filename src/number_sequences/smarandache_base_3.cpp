/*
* Creation date: 2026.05.19
* Created by Youcef Lemsafer
*/
#include "smarandache_base_3.hpp"

#include <string>
#include <cmath>

#include "lazy_init_array.hpp"

namespace cutrialdive {

    namespace {

        template <typename ValueT>
        ValueT compute_sm3_c_k(uint64_t k)
        {
            using IntT = number_sequence<num_seq_id::smarandache, 3>::value_type;
            // k is floor(log3(n)), 3^40 fits in an uint64_t
            constexpr auto MAX_PRECOMPUTED_CONSTANTS = 
#if CUTRIALDIVE_CPP_HAS_CONSTEXPR_LOGARITHM
                sizeof(number_sequence<num_seq_id::smarandache, 3>::index_type) * 8 * log(2) / log(3) + 1
#else
                40 + 1
#endif
            ;
            static auto firstValue = [](){
                if constexpr(std::is_same_v<std::string, ValueT>) {
                    return std::string{"3"};
                } else {
                    return ValueT{3};
                }
            };
            static lazy_init_array<ValueT, MAX_PRECOMPUTED_CONSTANTS> constants{firstValue(),
                [](ValueT const & predecessor, size_t i){
                // Remembering the following may ease the pain
                // c_i = (c_{i - 1}*3^(2*i*3^(i - 1)) - (3^i - 1)^2  - 3^i
                //       )/(3^i - 1)^2
                //        *(3^(i + 1) - 1)^2
                //       + 3^(2*i + 1) - 3^i + 1
                if constexpr(std::is_same_v<std::string, ValueT>) {
                    auto three_to_i = "3^" + std::to_string(i);
                    auto three_to_i_plus_1 = "3^" + std::to_string(i + 1);
                    auto three_to_two_i_plus_1 = "3^" + std::to_string(2*i + 1);
                    return "((" + predecessor + ")*3^" + to_string(2 * i * pow(IntT{3},i - 1)) + " - (" + three_to_i + " - 1)^2 - " + three_to_i
                            + ")/(" + three_to_i
                                + " - 1)^2*(" + three_to_i_plus_1 + " - 1)^2 + "
                                + three_to_two_i_plus_1 + " - " + three_to_i + " + 1";
                } else {
                    auto three_to_i = pow(IntT{3}, i);
                    auto three_to_i_minus_1_squared = pow(three_to_i - 1, 2);
                    return (predecessor * pow(three_to_i, 2 * pow(3, i - 1)) - three_to_i_minus_1_squared - three_to_i)
                            /three_to_i_minus_1_squared
                            *pow(three_to_i*3 - 1, 2)
                            + pow(three_to_i, 2)*3 - three_to_i + 1;
                }
            }};
            return constants[k];
        }

        constexpr uint64_t powersOf3[] = {1, 3, 9, 27, 81, 243, 729, 2187, 6561,
            19683, 59049, 177147, 531441, 1594323, 4782969, 14348907, 43046721,
            129140163, 387420489, 1162261467, 3486784401, 10460353203, 31381059609,
            94143178827, 282429536481, 847288609443, 2541865828329, 7625597484987,
            22876792454961, 68630377364883, 205891132094649, 617673396283947, 1853020188851841,
            5559060566555523, 16677181699666569, 50031545098999707, 150094635296999121,
            450283905890997363, 1350851717672992089, 4052555153018976267, 12157665459056928801u
        };

    } // anonymous namespace

    uint64_t log3(uint64_t n)
    {
        uint64_t res{};
        for(; n /= 3; ++res);
        return res;
    }

    char const* number_sequence<num_seq_id::smarandache, 3>::short_name()
    {
        return "Sm3";
    }

    void number_sequence<num_seq_id::smarandache, 3>::print_value(index_type n, std::ostream & out)
    {
        out << value(n);
    }


    typename number_sequence<num_seq_id::smarandache, 3>::value_type
    number_sequence<num_seq_id::smarandache, 3>::value(index_type n)
    {
        if(n < 2) {
            return n;
        }
        auto k = log3(n);
        auto c_k = compute_sm3_c_k<value_type>(k);
        // (C_k*3^((k + 1)*(n - 3^k + 1)) - (3^(k + 1) - 1)*n  - 3^(k + 1))/(3^(k + 1) - 1)^2
        auto three_to_k = powersOf3[k];
        return (c_k * pow(pow(value_type{3}, k + 1), n - three_to_k + 1) - (value_type{three_to_k} * 3 - 1)*n - value_type{three_to_k} * 3
                ) / pow(value_type{three_to_k} * 3 - 1, 2);
    }

    void number_sequence<num_seq_id::smarandache, 3>::print_expression(index_type n, std::ostream & out)
    {
        if(n < 2) {
            out << n;
            return;
        }
        auto k = log3(n);
        auto c_k = compute_sm3_c_k<std::string>(k);
        out << "((" << c_k << ")*3^" << (value_type{k + 1} *  value_type{n - powersOf3[k] + 1})
                    << " - (3^" << (k + 1) << " - 1)*" << n << " - 3^" << (k + 1)
                    << ")/(3^" << (k + 1) << " - 1)^2";
    }
}
