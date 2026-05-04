/*
* Created on 2026.05.03
* Copyright (c) Youcef Lemsafer
*/
#pragma once

#include <cstdint>

#include "common_defines.h"

namespace cutrialdive {

    template<typename ValueT>
    class device_span
    {
    public:
        device_span(ValueT * data, uint64_t size);
        CUTRIALDIVE_DEVICE ValueT& operator[](uint64_t i);
        CUTRIALDIVE_DEVICE ValueT const& operator[](uint64_t i) const;
    private:
        ValueT * data_;
        uint64_t size_;
    };
}

namespace cutrialdive {

    template <typename ValueT>
    inline
    device_span<ValueT>::device_span(ValueT * data, uint64_t size)
        : data_(data)
        , size_(size)
    {
    }

    template <typename ValueT>
    CUTRIALDIVE_DEVICE __forceinline__
    ValueT& device_span<ValueT>::operator[](uint64_t i)
    {
        return data_[i];
    }

    template <typename ValueT>
    CUTRIALDIVE_DEVICE __forceinline__
    ValueT const& device_span<ValueT>::operator[](uint64_t i) const
    {
        return data_[i];
    }

}
