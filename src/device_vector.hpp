/*
* MIT License
* Created on 2026.05.03
* Copyright (c) Youcef Lemsafer
* See LICENSE file for details.
*/
#pragma once

#include <cstdint>
#include "device_memory.hpp"
#include "device_span.hpp"

namespace cutrialdive {

    /// vector allocated on the device and managed on the host-side
    template <typename ValueT>
    class device_vector
    {
    public:
        device_vector() = default;
        device_vector(device_vector&&) = default;
        device_vector(device_vector const &) = delete;
        device_vector& operator=(device_vector const &) = delete;
        ~device_vector();

        /// Returns a device compatible view of the data
        device_span<ValueT> span() const;
        /// Returns pointer to data
        ValueT const * data() const;
        ValueT * data();
        /// Returns allocated size
        std::size_t size() const;
        /// Resize to @param newSize (does not reallocate if newSize <= size())
        void resize(std::size_t newSize);
        /// Resize to @param size and fill with zeroes
        void assign_zero(std::size_t size);

        std::vector<ValueT> to_host() const;

    private:
        ValueT * values_{};
        std::size_t size_{};
    };

}

namespace cutrialdive {

    template <typename ValueT>
    inline
    device_vector<ValueT>::~device_vector()
    {
        if(values_) {
            free_device_ptr(values_);
            values_ = nullptr;
        }
        size_ = 0;
    }

    template <typename ValueT>
    inline
    ValueT const * device_vector<ValueT>::data() const
    {
        return values_;
    }
    template <typename ValueT>
    inline
    ValueT * device_vector<ValueT>::data()
    {
        return values_;
    }

    template <typename ValueT>
    inline
    std::size_t device_vector<ValueT>::size() const
    {
        return size_;
    }

    template <typename ValueT>
    inline
    void device_vector<ValueT>::resize(std::size_t newSize)
    {
        if(newSize > size_) {
            free_device_ptr(values_);
            values_ = alloc_device_ptr<ValueT>(newSize);
        }
        size_ = newSize;
    }

    template <typename ValueT>
    inline
    void device_vector<ValueT>::assign_zero(std::size_t size)
    {
        resize(size);
        device_memset(values_, ValueT{0}, size_);
    }

    template <typename ValueT>
    inline
    device_span<ValueT> device_vector<ValueT>::span() const
    {
        return {values_, size_};
    }

    template <typename ValueT>
    inline
    std::vector<ValueT> device_vector<ValueT>::to_host() const
    {
        std::vector<ValueT> vec(size_, ValueT{0});
        copy_from_device(&vec[0], values_, size_);
        return vec;
    }

}
