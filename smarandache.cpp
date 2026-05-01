/*
 * Creation date: 2026.04.21
 * Created by Youcef Lemsafer
 */
#include "smarandache.hpp"

#include <algorithm>
#include <sstream>
#include <fstream>
#include <vector>
#include <memory>
#include <omp.h>

#include "hgint.hpp"
#include "siever.hpp"
#include "timer.hpp"


namespace cutrialdive
{
    ///
    std::ostream &output_smarandache(std::ostream &out, uint64_t n)
    {
        if(!n) {
            out << "0";
            return out;
        }
        for(decltype(n) i = 1; i <= n; ++i) {
            out << i;
        }
        return out;
    }

    std::ostream &output_smarandache_expression(std::ostream &out, uint64_t n)
    {
        if(n <= 1) {
            return out << n;
        }
        auto log10n = details::floor_log10(n);
        static const auto c3 = std::string{"1490840987654358*10^180-10990220"};
        static const auto c4 = std::string{"10^2701*("} + c3 + ")*10201-1109890222000";
        static const auto c5 = std::string{"10^36000*("} + c4 + ")*123454321-123443211358*33300^2";
        static const std::array<std::function<std::string(std::string)>, 5> sm = {
            [](auto x) { return std::string{"(10^("} + x + "+1)"
                                        + "-(9*" + x + "+10))/81"; }
          , [](auto x) { return std::string{"(10^(2*"} + x + "-18)*(123456789*99^2+991)"
                            "-(99*" + x + "+100))/9801"; }
          , [](auto x) { return std::string{"(("} + c3 + ")*10^(3*" + x + "-296)"
                            "-(120879*" + x + "+121000))/120758121"; }
          , [](auto x) { return std::string{"(("} + c4 + ")*10^(4*(" + x + "-999))"
                            "-12321*(9999*" + x + "+10000))/1231853592321"; }
          
          , [](auto x) { return std::string{"(("} + c5 + ")*10^(5*(" + x + "-9999))"
                                              "-12321*1234321*(99999*" + x + "+100000))"
                                            "/(12321*1234321*99999^2)"; }
        };
        if(log10n < sm.size()) {
            return out << sm[log10n](std::to_string(n));
        }
        throw std::runtime_error{"smarandache(n) for n >= 10^5 not supported yet"};
    }

    void read_previous_results(uint64_t n0, uint64_t n1, std::vector<result<uint64_t>> & previous_results)
    {
        auto n = n0;
        for(auto res : previous_results) {
            if(res.n != n) {
                std::ostringstream msg;
                msg << "Index mismatch between results and current range (from results: "
                    << res.n << ", from range: " << n << ")";
                throw std::runtime_error(msg.str());
            }
            ++n;
        }
    }

    template <typename T>
    void dump_results(std::vector<result<T>> const & results, std::ostream & out)
    {
        std::for_each(std::begin(results), std::end(results), [&out](auto const & res){
            out << res.n << " - ";
            char const * sep = "";
            std::for_each(std::begin(res.factors), std::end(res.factors), [&sep, &out](auto factor){
                out << sep << factor.value;
                if(factor.exponent > 1) {
                    out << "^" << factor.exponent;
                }
                sep = ", ";
            });
            out << "\n";
        });
    }

    uint64_t mu(uint64_t p)
    {
        return ~uint64_t{0} / p;
    }

    void do_trial_factor(uint64_t n0, uint64_t n1,
        uint64_t f0, uint64_t f1,
        std::vector<result<uint64_t>> & results
        )
    {
        constexpr auto max_num_of_primes = size_t{1} << 22;
        constexpr auto segment_len = uint64_t{1} << 26;

        std::vector<uint64_t> primes, inv_primes;
        primes.reserve(max_num_of_primes);
        inv_primes.reserve(max_num_of_primes);
        // residue of Sm(n0) modulo each prime
        std::vector<uint64_t> n0_residues;

        auto const sm_n0 = smarandache<HgInt>(n0);

        auto init_time = std::chrono::duration<double>{},
             propagation_time = init_time;

        for(auto fx = f0, fy = (std::min)(fx + segment_len, f1);
            fy <= f1;
            fx = fy, fy = (std::max)(fx + segment_len, f1)
        ) {
            primes.clear();
            inv_primes.clear();
            sieve(fx, fy, primes);
            if(n0_residues.size() < primes.size()) {
                n0_residues.reserve(primes.size());
                n0_residues.assign(primes.size(), uint64_t{});
            } else if(n0_residues.size() > primes.size()) {
                n0_residues.resize(primes.size());
            }
            inv_primes.assign(primes.size(), 0);

            // Compute Sm(n0) mod p and, if p != 2, the inverse of p mod 2^64.
            auto startTime = std::chrono::high_resolution_clock::now();
            #pragma omp parallel
            {
                HgInt local_sm_n0{sm_n0};
                size_t const i_max = primes.size();
                #pragma omp for
                for(size_t i = 0; i < i_max; ++i) {
                    auto p = primes[i];
                    n0_residues[i] = local_sm_n0.mod(p);
                    if(!n0_residues[i]) {
                        #pragma omp critical
                        results[0].factors.push_back({p, 1});
                    }
                    // We still push a value when p is 2 to avoid having the indice of
                    // p and mu(p) differ by one.
                    inv_primes[i] = p == 2 ? 0 : mu(p);
                }
            }
            auto now = std::chrono::high_resolution_clock::now();
            init_time += std::chrono::duration<double>(now - startTime);
            startTime = now;
            auto n = n0 + 1;
            for(auto i = 1; i < n1 - n0; ++i, ++n) {
                // 2 is so special that it dictates a special treatment.
                // Sieve results are always ordered so 2 can only be at the front.
                bool have_to_skip_front = false;
                if(!primes.empty() && primes.front() == 2) {
                    n0_residues[0] = n & 1;
                    if(!n0_residues[0]) {
                        results[i].factors.push_back({2, 1});
                    }
                    have_to_skip_front = true;
                }
                size_t j0 = have_to_skip_front ? 1 : 0, j_max = primes.size();
                #pragma omp parallel for
                for(size_t j = j0;
                    j < j_max;
                    ++j) {
                    auto a = __uint128_t(n0_residues[j]) * 100000 + n;
                    auto q = (a * inv_primes[j]) >> 64;
                    n0_residues[j] = a - q * primes[j];
                    if(n0_residues[j] >= primes[j]) {
                        n0_residues[j] -= primes[j];
                    }
                    if(!n0_residues[j]) {
                        #pragma omp critical
                        results[i].factors.push_back({primes[j], 1});
                    }
                }
                std::sort(std::begin(results[i].factors), std::end(results[i].factors),
                    [](auto const & x, auto const & y){ return x.value < y.value; });
            }
            now = std::chrono::high_resolution_clock::now();
            propagation_time += std::chrono::duration<double>(now - startTime);
        }
        std::cout << "Initialization time: " << init_time.count() << "s" << std::endl;
        std::cout << "Propagation time: " << propagation_time.count() << "s" << std::endl;
    }

    /// @brief Performs trial factoring according to given options
    /// @param opts trial factoring options
    void trial_factor(trial_factoring_options const & opts)
    {
        std::unique_ptr<std::ofstream> output_ptr;
        if(opts.output_path) {
            auto path = opts.output_path.value();
            if(std::filesystem::exists(path)) {
                std::string msg{"Output: cannot overwrite existing file "};
                msg += path.string();
                throw std::runtime_error(msg);
            }
            output_ptr = std::make_unique<std::ofstream>(path);
        }
        std::vector<result<uint64_t>> results;
        std::cout << "Trial factoring starts, results will be written to "
            << (output_ptr ? (std::string{"file `"} + opts.output_path.value().string() + "'") : "the console") << std::endl;
        time("Trial factoring took ", [&](){
            trial_factor(opts.n0, opts.n1, opts.f0, opts.f1, results);
        });
        dump_results(results, output_ptr ? *output_ptr.get() : std::cout);
    }


    /// @brief Perform the trial factoring of {Sm(i)} for i in [n0, n1[ using the primes
    /// in [f0, f1[.
    /// @param n0 lower bound of the range of indices to consider for the Smarandache numbers
    /// @param n1 upper bound of the range of indices to consider for the Smarandache numbers
    /// @param f0 lower bound of the set of primes
    /// @param f1 upper bound of the set of primes (f1 is excluded)
    /// @param previous_results previous results (can be empty e.g. first run)
    void trial_factor(uint64_t n0, uint64_t n1,
        uint64_t f0, uint64_t f1,
        std::vector<result<uint64_t>> & previous_results
    )
    {
        std::cout << "Maximum number of threads: " << omp_get_max_threads() << std::endl;
        if((n1 <= n0) || (f1 <= f0)) {
            return;
        }
        if(previous_results.empty()) {
            previous_results.reserve(n1 - n0);
            for(auto i = n0; i < n1; ++i) {
                previous_results.emplace_back(result<uint64_t>{i, {}});
            }
        } else {
            read_previous_results(n0, n1, previous_results);
        }
        do_trial_factor(n0, n1, f0, f1, previous_results);
    }


    void run_prp_test(uint64_t n, std::vector<factor<uint64_t>> & factors, bool haveToBoostFactors)
    {
        std::string const numberAsStr = std::string{"Sm("} + std::to_string(n) + ")";
        HgInt cofactor;
        {
            auto number = smarandache<HgInt>(n);
            timer t{(std::string{"Computing the cofactor of "}
                        + numberAsStr + " took ").c_str(), std::cout};
            std::cout << "Boosting factors: " << std::boolalpha << haveToBoostFactors << std::endl;
            cofactor = haveToBoostFactors ? compute_cofactor_boosted(number, factors)
                        : compute_cofactor_exact(number, factors);
        }
        {
            auto sizeInBase10 = cofactor.sizeInBase(10);
            timer t{(std::string{"The PRP test on the " + std::to_string(sizeInBase10) 
                + " digits cofactor of "} + numberAsStr
                + " took ").c_str(), std::cout};
            auto isPrp = is_prp(std::move(cofactor));
            std::cout << "Is PRP: " << (isPrp ? "YES" : "NO") << std::endl;
        }
    }
}
