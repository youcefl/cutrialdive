/*
* Created on 2026.04.21
* Copyright (c) Youcef Lemsafer
*/
#include "smarandache.hpp"

#include <string>
#include <sstream>
#include <limits>
#include <algorithm>

#include "lazy_init_array.hpp"


namespace cutrialdive {

    namespace {

        template <uint64_t base>
        std::vector<uint64_t> const & powers_of()
        {
            static auto const result = []() {
                std::vector<uint64_t> powersOfBase{1};
                for(auto e = base; ;) {
                    if(e > ~uint64_t{} / base) {
                        break;
                    }
                    powersOfBase.push_back(e);
                    e *= base;
                }
                return powersOfBase;
            }();
            return result;
        }

        template <uint64_t base, typename ValueT>
        ValueT compute_c_k(uint64_t k)
        {
            using IntT = typename smarandache<base>::value_type;
            using IndexT = typename smarandache<base>::index_type;
            constexpr auto indexBits = std::numeric_limits<IndexT>::digits;

            // Since k is floor(logb(n)), k < logb(2^indexBits) (where logb is the base b logarithm)
            // so the max value is indexBits - 1 (reached when base = 2), we could do better
            // when base > 2 (e.g. if base = 4 the max is actually indexBits / 2 - 1) but is it worth it?
            constexpr auto MAX_PRECOMPUTED_CONSTANTS = indexBits;
            static lazy_init_array<ValueT, MAX_PRECOMPUTED_CONSTANTS> constants{
                [](){
                    if constexpr(std::is_same_v<std::string, ValueT>) {
                        return to_string(base);
                    } else {
                        return ValueT{base};
                    }
                }(),
                [](ValueT const & predecessor, size_t i){

                // Remembering the following may ease the pain
                // C_i = (C_{i - 1}*b^((b - 1)*i*b^(i - 1))
                //       - (b^i - 1)^2
                //       - b^i
                //       )/(b^i - 1)^2
                //       *(b^(i + 1) - 1)^2 + b^(2*i + 1) - b^i + 1
                if constexpr(std::is_same_v<std::string, ValueT>) {
                    static const auto b = to_string(base);
                    auto b_to_i = (i != 1) ? b + "^" + std::to_string(i) : b;
                    auto b_to_i_plus_1 = b + "^" + std::to_string(i + 1);
                    auto b_to_two_i_plus_1 = b + "^" + std::to_string(2*i + 1);
                    auto pred_par = std::any_of(std::begin(predecessor), std::end(predecessor),
                        [](auto x){
                            return (x < '0') || (x > '9');
                        }) ? std::to_array<std::string>({"(", ")"})
                           : std::to_array<std::string>({"", ""});

                    return "("  + pred_par[0] + predecessor + pred_par[1] + "*" + b + "^" 
                                    + to_string((base - 1) * i * pow(IntT{base}, i - 1))
                                    + " - (" + b_to_i + " - 1)^2 - " + b_to_i
                            + ")/(" + b_to_i
                                + " - 1)^2*(" + b_to_i_plus_1 + " - 1)^2 + "
                                + b_to_two_i_plus_1 + " - " + b_to_i + " + 1";
                } else {
                    if(i >= 1 + powers_of<base>().size()) {
                        std::ostringstream ostr;
                        ostr << "While computing C_" << i << " for Smarandache base " 
                                << base << ", the value " << base << "^" << i 
                                << " which does not fit in 64 bits was requested";
                        throw std::runtime_error{ostr.str()};
                    }
                    constexpr auto b = base;
                    auto b_to_i = pow(IntT{b}, i);
                    auto b_to_i_minus_1_squared = pow(b_to_i - 1, 2);
                    return (predecessor * pow(pow(b_to_i, (b - 1)), powers_of<base>()[i - 1]) - b_to_i_minus_1_squared - b_to_i)
                            / b_to_i_minus_1_squared
                            * pow(b_to_i*b - 1, 2)
                            + pow(b_to_i, 2)*b - b_to_i + 1;
                }
            }};
            return constants[k];
        }

        template <uint64_t base>
        uint64_t logb(uint64_t n)
        {
            uint64_t res{};
            for(; n /= base; ++res);
            return res;
        }

    } // anonymous namespace

    template <uint64_t Base>
    char const* smarandache<Base>::short_name()
    {
        if constexpr (base == 10) {
            return "Sm";
        } else {
            static const auto name = std::string{"Sm"} + to_string(base);
            return name.c_str();
        }
    }

    template <uint64_t Base>
    std::ostream & smarandache<Base>::print_value(index_type n, std::ostream & out)
    {
        return out << value(n);
    }

    template <uint64_t Base>
    typename smarandache<Base>::value_type
    smarandache<Base>::value(index_type n)
    {
        if(n < 2) {
            return n;
        }
        auto k = logb<base>(n);
        auto c_k = compute_c_k<base, value_type>(k);
        // b being the base, recall that
        // Smb(n) = (C_k*b^((k + 1)*(n - b^k + 1))
        //          - (b^(k + 1) - 1)*n
        //          - b^(k + 1)
        //          )/(b^(k + 1) - 1)^2
        if constexpr(std::has_single_bit(base)) {
            // When base is a power of 2
            constexpr auto shift_multiplier = floor_log2_base;
            return ((c_k << (shift_multiplier*(k + 1)*(n - (1 << (shift_multiplier * k)) + 1)))
                    - (((1 << (shift_multiplier * k)) << shift_multiplier) - 1)*n
                    - ((1 << (shift_multiplier * k)) << shift_multiplier)
                    )/pow(((1 << (shift_multiplier * k)) << shift_multiplier) - 1, 2);
        } else {
            static const auto b = value_type{base};
            auto b_to_k_plus_one = pow(b, k + 1);
            return ((c_k*pow(b_to_k_plus_one, n - powers_of<base>()[k] + 1))
                    - (b_to_k_plus_one - 1)*n
                    - b_to_k_plus_one
                   )/pow(b_to_k_plus_one - 1, 2);
        }
    }

    template <uint64_t Base>
    std::ostream & smarandache<Base>::print_expression(index_type n, std::ostream & out)
    {
        if(n < 2) {
            return out << n;
        }
        auto k = logb<base>(n);
        auto c_k = compute_c_k<base, std::string>(k);
        if(std::any_of(std::begin(c_k), std::end(c_k), [](auto x){
            return (x < '0') || (x > '9');
        })) {
            c_k = "(" + c_k + ")";
        }
        auto base_to_k_plus_1 = (k == 0)
                                ? to_string(base) 
                                : to_string(base) + "^" + to_string(k + 1);
        out << "(" << c_k << "*" << base << "^" << (value_type{k + 1} * (n - powers_of<base>()[k] + 1))
            << " - (" << base_to_k_plus_1 << " - 1)*" << n << " - " << base_to_k_plus_1
            << ")/(" << base_to_k_plus_1 << " - 1)^2";
        return out;
    }

    template struct smarandache<2>;
    template struct smarandache<3>;
    template struct smarandache<4>;
    template struct smarandache<5>;
    template struct smarandache<6>;
    template struct smarandache<7>;
    template struct smarandache<8>;
    template struct smarandache<9>;
    template struct smarandache<10>;
    template struct smarandache<11>;
    template struct smarandache<12>;
    template struct smarandache<13>;
    template struct smarandache<14>;
    template struct smarandache<15>;
    template struct smarandache<16>;
    template struct smarandache<17>;
    template struct smarandache<18>;
    template struct smarandache<19>;
    template struct smarandache<20>;
    template struct smarandache<21>;
    template struct smarandache<22>;
    template struct smarandache<23>;
    template struct smarandache<24>;
    template struct smarandache<25>;
    template struct smarandache<26>;
    template struct smarandache<27>;
    template struct smarandache<28>;
    template struct smarandache<29>;
    template struct smarandache<30>;
    template struct smarandache<31>;
    template struct smarandache<32>;
    template struct smarandache<33>;
    template struct smarandache<34>;
    template struct smarandache<35>;
    template struct smarandache<36>;
    template struct smarandache<37>;
    template struct smarandache<38>;
    template struct smarandache<39>;
    template struct smarandache<40>;
    template struct smarandache<41>;
    template struct smarandache<42>;
    template struct smarandache<43>;
    template struct smarandache<44>;
    template struct smarandache<45>;
    template struct smarandache<46>;
    template struct smarandache<47>;
    template struct smarandache<48>;
    template struct smarandache<49>;
    template struct smarandache<50>;
    template struct smarandache<51>;
    template struct smarandache<52>;
    template struct smarandache<53>;
    template struct smarandache<54>;
    template struct smarandache<55>;
    template struct smarandache<56>;
    template struct smarandache<57>;
    template struct smarandache<58>;
    template struct smarandache<59>;
    template struct smarandache<60>;
    template struct smarandache<61>;
    template struct smarandache<62>;
    template struct smarandache<63>;
    template struct smarandache<64>;
}

