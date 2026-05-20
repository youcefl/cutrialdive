/*
* Created on 2026.05.15
* Copyright (c) Youcef Lemsafer
*/
#pragma once

#include <string>

#include "num_seq_id.hpp"

namespace cutrialdive {

    /// Number sequence specification
    struct num_seq_spec
    {
        num_seq_id num_seq_id_value;
        std::string num_seq_params;
    };
}
