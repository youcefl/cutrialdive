/*
* Created on 2026.06.05
* Copyright (c) Youcef Lemsafer
*/
#include "progress.hpp"
#include "modular_arithmetic_detail.hpp"

#include <iomanip>
#include <string>


namespace {

    void print_duration(std::ostream & out, double seconds)
    {
        if(seconds > 365*86400) {
            out << std::setw(8)
                << std::fixed
                << std::setprecision(2)
                << seconds / (365 * 86400) << " years";
            return;
        }
        if (seconds > 86400) {
            out << std::setw(3)
                << std::fixed
                << int32_t(seconds / 86400) << "d ";
            seconds = seconds - 86400 * int32_t(seconds / 86400);
        }
        auto h = int32_t(seconds / 3600);
        auto m = int32_t((seconds - 3600*h) / 60);
        auto s = int32_t(seconds - 3600*h - 60*m);
        out << std::setw(2) << std::setfill('0') << h
            << ":" << std::setw(2) << std::setfill('0') << m
            << ":" << std::setw(2) << std::setfill('0') << s
            ;
        out << std::setfill(' ');
    }

    void print_eta(std::ostream & out, double percentDone, double elapsedTimeSeconds)
    {
        auto eta = (100 - percentDone)/percentDone * elapsedTimeSeconds;
        print_duration(out, eta);
    }

}

namespace cutrialdive {

    progress::progress(uint64_t start, uint64_t target, std::chrono::milliseconds period, std::ostream & out)
        : start_(start)
        , target_(target)
        , period_(period)
        , out_(out)
        , ended_(false)
        , target_digits_(digits_in_base<10>(target_))
        , start_time_(std::chrono::steady_clock::now())
        , last_update_time_(start_time_)
    {
    }

    progress::~progress()
    {
        if(!ended_) {
            out_ << std::endl;
        }
    }

    void progress::update(uint64_t current, uint64_t lastPrime)
    {
        constexpr auto inFactor = 10'000'000'000;
        constexpr auto outFactor = 100'000'000.0;

        auto now = std::chrono::steady_clock::now();
        if((std::chrono::duration_cast<std::chrono::milliseconds>(now - last_update_time_) < period_)
          && (current != target_)) {
            return;
        }
        auto percentDone = double(((__uint128_t(current) - start_) * inFactor) / (target_ - start_)) / outFactor;
        auto elapsedTime = std::chrono::duration<double>(now - start_time_).count();

        static const std::string spaces(128, ' ');

        out_ << "\r" << spaces
            << "\r["
            << std::setw(6)
            << std::fixed
            << std::setprecision(2)
            << percentDone << "% done]"
            << "["
            << std::setw(target_digits_)
            << lastPrime
            << "][Elapsed ";
        print_duration(out_, elapsedTime);
        out_ << "]";
        if(current != target_) {
            out_ << "[ETA ";
            print_eta(out_, percentDone, elapsedTime);
            out_ << "]";
        }
        out_ << std::flush;
        last_update_time_ = std::chrono::steady_clock::now();
    }

    void progress::end()
    {
        if(ended_) {
            return;
        }
        out_ << std::endl;
        ended_ = true;
    }
}
