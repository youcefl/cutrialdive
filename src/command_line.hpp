/*
* Created on 2026.05.04
* Copyright (c) Youcef Lemsafer
*/
#pragma once

#include <iostream>
#include <optional>

#include "trial_factoring_options.hpp"
#include "factor.hpp"
#include "mode.hpp"

namespace cutrialdive {

    void print_usage(std::ostream & out);
    void print_help(std::ostream & out);
    void print_version(std::ostream & out);

    struct command_line_options
    {
        mode_flag mode;
        bool wants_value;
        bool wants_expression;
        uint64_t n;
        bool wants_single_prp;
        bool wants_boosted_factors;
        std::vector<factor<uint64_t, uint32_t>> factors;
        std::optional<trial_factoring_options> tf_options;
    };
    
    class command_line_parser
    {
    public:
        command_line_parser(int argc, char** argv);
        bool is_help_requested() const;
        bool is_usage_requested() const;
        bool is_version_requested() const;
        std::optional<command_line_options> get_options() const;
    private:
        bool wants_help_;
        bool wants_usage_;
        bool wants_version_;
        std::optional<command_line_options> options_;
    };

    class command_line_parse_exception : public std::runtime_error
    {
    private:
        command_line_parse_exception(std::string const & msg);
        command_line_parse_exception(char const* msg);
        friend class command_line_parser;
    };
}


namespace cutrialdive {

    inline
    bool command_line_parser::is_help_requested() const
    {
        return wants_help_;
    }
    inline
    bool command_line_parser::is_usage_requested() const
    {
        return wants_usage_;
    }
    inline
    bool command_line_parser::is_version_requested() const
    {
        return wants_version_;
    }
    inline
    std::optional<command_line_options> command_line_parser::get_options() const
    {
        return options_;
    }

}