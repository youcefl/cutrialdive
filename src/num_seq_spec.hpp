/*
* MIT License
* Created on 2026.05.15
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#pragma once

#include <string>

#include "num_seq_id.hpp"

namespace cutrialdive {

    /// Number sequence specification
    struct num_seq_spec
    {
        num_seq_id seq_id;
        std::string seq_params;

        bool operator==(num_seq_spec const&) const = default;        
    };

}
