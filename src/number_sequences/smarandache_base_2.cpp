/*
 * Creation date: 2026.05.15
 * Created by Youcef Lemsafer
 */
#include "smarandache_base_2.hpp"

#include "lazy_init_array.hpp"
#include <array>
#include <sstream>
#include <functional>

// See doc at the end of the file on how the formula is derived

namespace cutrialdive {

    namespace {

        std::string mersenne_as_string(uint64_t n)
        {
            using IntT = number_sequence<num_seq_id::smarandache, 2>::value_type;
            return to_string((IntT{1} << n) - 1);
        }

        template <typename ValueT>
        ValueT compute_sm2_c_k(uint64_t k)
        {
            using IntT = number_sequence<num_seq_id::smarandache, 2>::value_type;
            // When Sm2(n) is requested this function is called with
            // k = floor(log2(n)) - 1 = std::bit_width(n) - 2
            // n being a 64-bit value
            constexpr auto MAX_PRECOMPUTED_CONSTANTS = 64 - 2 + 1;
            static auto firstValue = [](){
                if constexpr(std::is_same_v<std::string, ValueT>) {
                    return std::string{"16"};
                } else {
                    return ValueT{16};
                }
            };
            static lazy_init_array<ValueT, MAX_PRECOMPUTED_CONSTANTS> constants{firstValue(),
                [](ValueT const & predecessor, size_t i){
                if constexpr(std::is_same_v<std::string, ValueT>) {
                    auto mersenne_i_plus_1 = mersenne_as_string(i + 1);
                    auto mersenne_i_plus_2 = mersenne_as_string(i + 2);
                    auto delta = (IntT{1} << (2*i + 3)) - (IntT{1} << (i + 1)) + 1;
                    return mersenne_i_plus_2 + "^2*(("
                            + predecessor + ")*2^" + to_string(IntT{i + 1} << i)
                            + "-" + to_string(pow((IntT{1} << (i + 1)) - 1, 2) + (IntT{1} << (i + 1)))
                            + ")/" + mersenne_i_plus_1 + "^2"
                            "+" + to_string(delta);
                } else {
                    // @todo: the shift overflows if i >= 59, if that happens it means
                    // someone requested Sm2(n) for an n >= 2^60
                    auto successor = predecessor << ((i + 1) * (1ull << i));
                    auto mersenne_i_plus_1 = (ValueT{1} << (i + 1)) - 1;
                    successor = pow((mersenne_i_plus_1 << 1) + 1, 2) 
                            * (successor - pow(mersenne_i_plus_1, 2) - (ValueT{1} << (i + 1)))
                            / pow(mersenne_i_plus_1, 2)
                            + (ValueT{1} << (2 * i + 3)) - mersenne_i_plus_1;
                    return std::move(successor);
                }

            }};
            return constants[k];
        }
    } // anonymous namespace


    char const* number_sequence<num_seq_id::smarandache, 2>::short_name()
    {
        return "Sm2";
    }

    void number_sequence<num_seq_id::smarandache, 2>::print_value(index_type n, std::ostream & out)
    {
        out << value(n);
    }

    typename number_sequence<num_seq_id::smarandache, 2>::value_type
    number_sequence<num_seq_id::smarandache, 2>::value(index_type n)
    {
        if(n < 2) {
            return value_type{(n <= 0) ? 0u : 1u};
        }
        auto const k = std::bit_width(n) - 2;
        auto c_k = compute_sm2_c_k<number_sequence<num_seq_id::smarandache, 2>::value_type>(k);

        auto c_k_shift = (k + 2) * (n - (1ull << (k + 1)) + 1);
        auto mersenne_k_plus_2 = (value_type{1} << (k + 2)) - 1;
        return ((c_k << c_k_shift) - mersenne_k_plus_2 * n - mersenne_k_plus_2 - 1)
                / (mersenne_k_plus_2 * mersenne_k_plus_2);

    }

    void number_sequence<num_seq_id::smarandache, 2>::print_expression(index_type n, std::ostream & out)
    {
        if(n < 2) {
            out << ((n <= 0) ? "0" : "1");
            return;
        }
        auto const k = std::bit_width(n) - 2;
        auto c_k = compute_sm2_c_k<std::string>(k);

        auto c_k_shift = value_type{k + 2}*(value_type{n} - (value_type{1} << (k + 1)) + 1);
        auto mersenne_k_plus_2 = (value_type{1} << (k + 2)) - 1;
        out << "((" << c_k << ")*2^" << c_k_shift
            << "-" << (mersenne_k_plus_2 * n + mersenne_k_plus_2 + 1) << ")/"
            << mersenne_k_plus_2 << "^2";
    }

}


    /*
    Definition: Sm2(n) is Smarandache base 2 i.e. concatenation of the first n integers written in base 2
    Sm2(0) = 0
    Sm2(1) = 1
    Sm2(2) = 110_2 = 6
    Sm2(3) = 11011_2 = 27
    Sm2(4) = 11011100_2 = 220

    Let c0 = 1 = Sm2(1), for 2 <= n < 4 we have
        Sm2(n) = c0*2^(2*(n - 1)) + sum(k=2, n, k*2^(2*(n - k)))
            = c0*2^(2*(n - 1)) + (7*2^(2*(n - 1)) - 3*n - 4)/9
            = ((9*c0 + 7)*2^(2*(n - 1)) - 3*n - 4)/9
    Let c1 = Sm2(3) = ((9*c0 + 7)*2^4 - 13)/9, for 4 <= n < 8 we have
        Sm2(n) = c1*2^(3*(n-3)) + sum(k=4, n, k*2^(3*(n-k)))
            = c1*2^(3*(n-3)) + (29*2^(3*n - 9) - 7*n - 8)/49
            = ((49*c1 + 29)*2^(3*(n - 3)) - 7*n - 8)/49
    Let c2 = Sm2(7) = ((49*c1 + 29)*2^12 - 57)/49, for 8 <= n < 16 we have
        Sm2(n) = c2*2^(4*(n-7)) + sum(k=8, n, k*2^(4*(n-k)))
            = c2*2^(4*(n-7)) + (121*2^(4*(n - 7)) - 15*n - 16)/225
            = ((225*c2 + 121)*2^(4*(n - 7)) - 15*n - 16)/225
    Let c3 = Sm2(15) = ((225*c2 + 121)*2^32 - 241)/225, then for 16 <= n < 32 we have
        Sm2(n) = c3*2^(5*(n-15)) + sum(k=16, n, k*2^(5*(n-k)))
            = c3*2^(5*(n-15)) + (497*2^(5*(n - 15)) - 31*n - 32)/961
            = ((961*c3 + 497)*2^(5*(n - 15)) - 31*n - 32)/961
    Let c4 = Sm2(31) = ((961*c3 + 497)*2^80 - 993)/961, then for 32 <= n < 64 we have
        Sm2(n) = c4*2^(6*(n-31)) + sum(k=32, n, k*2^(6*(n-k)))
            = c4*2^(6*(n-31)) + (2017*2^(6*(n - 31)) - 63*n - 64)/3969
            = ((3969*c4 + 2017)*2^(6*(n - 31)) - 63*n - 64)/3969
    Let c5 = Sm2(63) = ((3969*c4 + 2017)*2^192 - 4033)/3969, then for 64 <= n < 128 we have
        Sm2(n) = c5*2^(7*(n - 63)) + sum(k=64, n, k*2^(7*(n - k)))
            = c5*2^(7*(n - 63)) + (8129*2^(7*(n - 63)) - 127*n - 128)/16129
            = ((c5*16129 + 8129)*2^(7*(n - 63)) - 127*n - 128)/16129
    Let c6 = Sm2(127) = ((c5*16129 + 8129)*2^448 - 16257)/16129, then for 128 <= n < 256 we have
        Sm2(n) = c6*2^(8*(n - 127)) + sum(k=128, n, k*2^(8*(n - k)))
               = c6*2^(8*(n - 127)) + (32641*2^(8*(n - 127)) - 255*n - 256)/65025
               = ((65025*c6 + 32641)*2^(8*(n - 127)) - 255*n - 256)/65025
    Let c7 = Sm2(255) = ((65025*c6 + 32641)*2^1024 - 65281)/65025, then for 256 <= n < 512 we have
        Sm2(n) = c7*2^(9*(n - 255)) + sum(k=256, n, k*2^(9*(n - k)))
               = c7*2^(9*(n - 255)) + (130817*2^(9*(n - 255)) - 511*n - 512)/261121
               = ((261121*c7 + 130817)*2^(9*(n - 255)) - 511*n - 512)/261121
    Let c8 = Sm2(511) = ((261121*c7 + 130817)*2^2304 - 261633)/261121, then for 512 <= n < 1024 we have
        Sm2(n) = c8*2^(10*(n - 511)) + sum(k=512, n, k*2^(10*(n - k)))
               = c8*2^(10*(n - 511)) + (523777*2^(10*(n - 511)) - 1023*n - 1024)/1046529
               = ((1046529*c8 + 523777)*2^(10*(n - 511)) - 1023*n - 1024)/1046529
    Let c9 = Sm2(1023) = ((1046529*c8 + 523777)*2^5120 - 1047553)/1046529, then for 1024 <= n < 2048 we have
        Sm2(n) = c9*2^(11*(n - 1023)) + sum(k=1024, n, k*2^(11*(n - k)))
               = c9*2^(11*(n - 1023)) + (2096129*2^(11*(n - 1023)) - 2047*n - 2048)/4190209
               = ((4190209*c9 + 2096129)*2^(11*(n - 1023)) - 2047*n - 2048)/4190209
    Let c10 = Sm2(2047) = ((4190209*c9 + 2096129)*2^11264 - 4192257)/4190209, then for 2048 <= n < 4096 we have
        Sm2(n) = c10*2^(12*(n - 2047)) + sum(k=2048, n, k*2^(12*(n - k)))
               = c10*2^(12*(n - 2047)) + (8386561*2^(12*(n - 2047)) - 4095*n - 4096)/4095^2
               = ((4095^2*c10 + 8386561)*2^(12*(n - 2047)) - 4095*n - 4096)/4095^2
    Let c11 = Sm2(4095) = ((4095^2*c10 + 8386561)*2^24576 - 16773121)/4095^2, then for 4096 <= n < 8192 we have
        Sm2(n) = c11*2^(13*(n - 4095)) + sum(k=4096, n, k*2^(13*(n - k)))
               = c11*2^(13*(n - 4095)) + (33550337*2^(13*(n - 4095)) - 8191*n - 8192)/8191^2
               = ((8191^2*c11 + 33550337)*2^(13*(n - 4095)) - 8191*n - 8192)/8191^2
    Let c12 = Sm2(8191) = ((8191^2*c11 + 33550337)*2^53248 - 67100673)/8191^2, then for 8192 <= n < 16384 we have
        Sm2(n) = c12*2^(14*(n - 8191)) + sum(k=8192, n, k*2^(14*(n - k)))
               = c12*2^(14*(n - 8191)) + (134209537*2^(14*(n - 8191)) - 16383*n - 16384)/16383^2
               = ((16383^2*c12 + 134209537)*2^(14*(n - 8191)) - 16383*n - 16384)/16383^2


    We have c0 = 16,
            c1 = 49*(c0*2^4 - 13)/9 + 29,
            c2 = 225*(c1*2^12 - 57)/49 + 121,
            c3 = 961*(c2*2^32 - 241)/225 + 497,
            c4 = 3969*(c3*2^80 - 993)/961 + 2017,
            c5 = 16129*(c4*2^192 - 4033)/3969 + 8129,
            c6 = 65025*(c5*2^448 - 16257)/16129 + 32641,
            c7 = 261121*(c6*2^1024 - 65281)/65025 + 130817,
            c8 = 1046529*(c7*2^2304 - 261633)/261121 + 523777,
            c9 = 4190209*(c8*2^5120 - 1047553)/1046529 + 2096129,
            c10 = 4095*4095*(c9*2^11264 - 4192257)/4190209 + 8386561,
            c11 = 8191*8191*(c10*2^24576 - 16773121)/(4095*4095) + 33550337,
            c12 = 16383*16383*(c11*2^53248 - 67100673)/(8191*8191) + 134209537

    So we can define c_k this way: c_0 = 16. For k >= 0
    c_{k+1} = (2^(k + 3) - 1)^2*(c_{k}*2^((k + 2)*2^(k + 1)) - (2^(k + 2) - 1)^2 - 2^(k + 2))/(2^(k + 2) - 1)^2
                    + 2^(2*k + 5) - 2^(k + 2) + 1

    We can then use PARI/GP to generate initialization code for the {c_i}:
      u(k) = 2^(k + 3) - 1
      v(k) = (k + 2)*2^(k + 1)
      w(k) = (2^(k + 2) - 1)^2 + 2^(k + 2)
      x(k) = 2^(k + 2) - 1
      y(k) = 2^(2*k + 5) - 2^(k + 2) + 1
      for(k = 0, 18, print("c", k+1, " = ", u(k), "ull*", u(k), "*((c", k,
            " << ", v(k), ") - ", w(k), "ull)/(", x(k), "ull*", x(k), ") + ", y(k), "ull,"))
    */

        // auto compute_sm2_c_i(size_t i)
        // {
        //     using Value = number_sequence<num_seq_id::smarandache, 2>::value_type;
        //     static Value c0{16},
        //         c1 = 7ull*7*((c0 << 4) - 13ull)/(3ull*3) + 29ull,
        //         c2 = 15ull*15*((c1 << 12) - 57ull)/(7ull*7) + 121ull,
        //         c3 = 31ull*31*((c2 << 32) - 241ull)/(15ull*15) + 497ull,
        //         c4 = 63ull*63*((c3 << 80) - 993ull)/(31ull*31) + 2017ull,
        //         c5 = 127ull*127*((c4 << 192) - 4033ull)/(63ull*63) + 8129ull,
        //         c6 = 255ull*255*((c5 << 448) - 16257ull)/(127ull*127) + 32641ull,
        //         c7 = 511ull*511*((c6 << 1024) - 65281ull)/(255ull*255) + 130817ull,
        //         c8 = 1023ull*1023*((c7 << 2304) - 261633ull)/(511ull*511) + 523777ull,
        //         c9 = 2047ull*2047*((c8 << 5120) - 1047553ull)/(1023ull*1023) + 2096129ull,
        //         c10 = 4095ull*4095*((c9 << 11264) - 4192257ull)/(2047ull*2047) + 8386561ull,
        //         c11 = 8191ull*8191*((c10 << 24576) - 16773121ull)/(4095ull*4095) + 33550337ull,
        //         c12 = 16383ull*16383*((c11 << 53248) - 67100673ull)/(8191ull*8191) + 134209537ull,
        //         c13 = 32767ull*32767*((c12 << 114688) - 268419073ull)/(16383ull*16383) + 536854529ull,
        //         c14 = 65535ull*65535*((c13 << 245760) - 1073709057ull)/(32767ull*32767) + 2147450881ull,
        //         c15 = 131071ull*131071*((c14 << 524288) - 4294901761ull)/(65535ull*65535) + 8589869057ull,
        //         c16 = 262143ull*262143*((c15 << 1114112) - 17179738113ull)/(131071ull*131071) + 34359607297ull,
        //         c17 = 524287ull*524287*((c16 << 2359296) - 68719214593ull)/(262143ull*262143) + 137438691329ull,
        //         c18 = 1048575ull*1048575*((c17 << 4980736) - 274877382657ull)/(524287ull*524287) + 549755289601ull,
        //         c19 = 2097151ull*2097151*((c18 << 10485760) - 1099510579201ull)/(1048575ull*1048575) + 2199022206977ull;
            

        // switch(idx) {
        //     case 0: return ((c0 << (2*(n - 1))) - 3*n - 4)/9;
        //     case 1: return ((c1 << (3*(n - 3))) - 7*n - 8)/49;
        //     case 2: return ((c2 << (4*(n - 7))) - 15*n - 16)/225;
        //     case 3: return ((c3 << (5*(n - 15))) - 31*n - 32)/961;
        //     case 4: return ((c4 << (6*(n - 31))) - 63*n - 64)/3969;
        //     case 5: return ((c5 << (7*(n - 63))) - 127*n - 128)/16129;
        //     case 6: return ((c6 << (8*(n - 127))) - 255*n - 256)/65025;
        //     case 7: return ((c7 << (9*(n - 255))) - 511*n - 512)/261121;
        //     case 8: return ((c8 << (10*(n - 511))) - 1023*n - 1024)/1046529;
        //     case 9: return ((c9 << (11*(n - 1023))) - 2047*n - 2048)/4190209;
        //     case 10: return ((c10 << (12*(n - 2047))) - 4095*n - 4096)/(4095*4095);
        //     case 11: return ((c11 << (13*(n - 4095))) - 8191*n - 8192)/(8191*8191);
        //     case 12: return ((c12 << (14*(n - 8191))) - 16383*n - 16384)/(16383*16383);
        //     default: {
        //         std::ostringstream ostr;
        //         ostr << "Cannot compute Sm2(n) when n >= 2^" << maxPowOf2 << " = " << (1ull << maxPowOf2);
        //         throw std::runtime_error{ostr.str()};
        //     }
        // }

    // Keep this somewhere:
    64 <= n < 128: ((127^2*(2329049770201871531016580615632863189041216*2^192-4033)/63^2+8129)*2^(7*(n-63))-127*n-128)/127^2
    128 <= n < 256:
    256 <= n < 512: ((511^2*((255^2*((127^2*(2329049770201871531016580615632863189041216*2^192-4033)/63^2+8129)*2^448-16257)/127^2+32641)*2^1024-65281)/255^2+130817)*2^(9*(n - 255))-511*n-512)/511^2
    512 <= n < 1024: ((1023^2*((511^2*((255^2*((127^2*(2329049770201871531016580615632863189041216*2^192-4033)/63^2+8129)*2^448-16257)/127^2+32641)*2^1024-65281)/255^2+130817)*2^2304-261633)/511^2+523777)*2^(10*(n - 511))-1023*n-1024)/1023^2

