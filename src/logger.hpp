/*
* MIT License
* Created on 2026.06.18
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#pragma once

#include <ostream>
#include <memory>
#include "verbosity.hpp"

#define CTDLOG_PASTE(u, v) u##v
#define CTDLOG_LOOP_VAR_NAME(l) CTDLOG_PASTE(cutrialdiveLogHaveToLoop, l)

#define CTDLOG_LOOP(line, v, f) \
    for(auto CTDLOG_LOOP_VAR_NAME(line) = true; CTDLOG_LOOP_VAR_NAME(line); CTDLOG_LOOP_VAR_NAME(line) = false) \
        if(log().get_verbosity() >= verbosity::v) \
            log().f()

#define CTDLOG_XXX2(v, f)  CTDLOG_LOOP(__LINE__, v, f)
#define CTDLOG_XXX(v) CTDLOG_XXX2(v, v)

#define CTDLOG_INFO()      CTDLOG_XXX2(normal, info)
#define CTDLOG_VERBOSE()   CTDLOG_XXX(verbose)
#define CTDLOG_DEBUG()     CTDLOG_XXX(debug)
#define CTDLOG_DIAG()      CTDLOG_XXX(diagnostic)

namespace cutrialdive {

    class logger
    {
    public:
        struct stream
        {
            std::ostream* output;
            
            template<typename T>
            stream& operator<<(T&& val);
            stream& operator<<(std::ostream& (*manip)(std::ostream&));
        };
    
        stream info();
        stream verbose();
        stream debug();
        stream diagnostic();

        verbosity get_verbosity() const;
        void set_verbosity(verbosity verbosityLevel);
        void set_stream(std::shared_ptr<std::ostream> out);
        void set_stream(std::ostream & out);


    private:
        logger(verbosity maxLevel, std::ostream & out);
        verbosity max_level_;
        std::shared_ptr<std::ostream> out_;

        friend logger & log();
    };

    logger & log();
}


namespace cutrialdive {

    inline
    logger::logger(verbosity maxLevel, std::ostream & out)
        : max_level_(maxLevel)
        , out_(std::shared_ptr<std::ostream>{&out, [](std::ostream*){}})
    {}

    inline
    void logger::set_verbosity(verbosity verbosityLevel)
    {
        max_level_ = verbosityLevel;
    }

    inline
    verbosity logger::get_verbosity() const
    {
        return max_level_;
    }

    inline
    void logger::set_stream(std::shared_ptr<std::ostream> out)
    {
        out_ = out;
    }

    inline
    void logger::set_stream(std::ostream & out)
    {
        out_ = std::shared_ptr<std::ostream>{&out, [](std::ostream*){}};
    }

    inline
    logger::stream logger::info()
    {
        return stream{max_level_ >= verbosity::normal ? out_.get() : nullptr};
    }

    inline
    logger::stream logger::verbose()
    {
        return stream{max_level_ >= verbosity::verbose ? out_.get() : nullptr};
    }

    inline
    logger::stream logger::debug()
    {
        return stream{max_level_ >= verbosity::debug ? out_.get() : nullptr};
    }

    inline
    logger::stream logger::diagnostic()
    {
        return stream{max_level_ >= verbosity::diagnostic ? out_.get() : nullptr};
    }

    template<typename T>
    inline
    logger::stream& logger::stream::operator<<(T&& val)
    {
        if(output) {
            *output << std::forward<T>(val);
        }
        return *this;
    }

    inline
    logger::stream& logger::stream::operator<<(std::ostream& (*manip)(std::ostream&))
    {
        if(output) {
            manip(*output);
        }
        return *this;
    }    
}
