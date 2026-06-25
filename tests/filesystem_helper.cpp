/*
* MIT License
* Created on 2026.06.09
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#ifdef _WIN32
#  include <windows.h>
#else
#  include <unistd.h>
#endif

#include <cstdint>
#include <array>
#include <random>
#include <string>
#include <iomanip>
#include <sstream>


#include "filesystem_helper.hpp"

namespace fs = std::filesystem;

namespace {

    inline uint64_t pid()
    {
#ifdef _WIN32
        return static_cast<uint64_t>(::GetCurrentProcessId());
#else
        return static_cast<uint64_t>(::getpid());
#endif
    }

    std::string timestamp_centis()
    {
        using namespace std::chrono;

        auto now = system_clock::now();
        auto secs = time_point_cast<seconds>(now);
        auto ms   = duration_cast<milliseconds>(now - secs);

        std::time_t t = system_clock::to_time_t(secs);

        std::tm tm{};
#ifdef _WIN32
        gmtime_s(&tm, &t);
#else
        gmtime_r(&t, &tm);
#endif

        std::ostringstream oss;

        oss << std::put_time(&tm, "%Y%m%dT%H%M%S")
            << std::setw(2) << std::setfill('0')
            << (ms.count() / 10);   // centiseconds

        return oss.str();
    }
 
    std::string random_12_bytes_hex()
    {
        std::array<uint8_t, 12> bytes;

        std::random_device rd;
        for (auto &b : bytes)
            b = static_cast<uint8_t>(rd());

        static constexpr char hex[] = "0123456789abcdef";

        std::string out;
        out.reserve(24);

        for (auto b : bytes)
        {
            out.push_back(hex[b >> 4]);
            out.push_back(hex[b & 0x0F]);
        }

        return out;
    }    
}

namespace cutrialdive::tests {



    fs::path get_temp_file_path()
    {
        // <tmppath>/<pid>_<16_random_bytes>
        
        return fs::temp_directory_path() 
               / (timestamp_centis() + "_" + std::to_string(pid()) + "_" + random_12_bytes_hex());
    }
}

