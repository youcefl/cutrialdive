/*
* Created on 2026.05.03
* Copyright (c) Youcef Lemsafer
*/
#pragma once

 #include <cstdint>

namespace cutrialdive {

    template <typename PtrT>
    PtrT* alloc_device_ptr(uint64_t size);

    template <typename PtrT>
    void free_device_ptr(PtrT *& ptr);

    template <typename ValueT>
    void device_memset(ValueT * dst, ValueT value, std::size_t size);

    template <typename ValueT>
    void copy_from_device(ValueT * dst, ValueT const * src, std::size_t size);

    template <typename ValueT>
    void copy_to_device(ValueT * dst, ValueT const * src, std::size_t size);
}

namespace cutrialdive {

    template <typename PtrT>
    inline
    void free_device_ptr(PtrT *& ptr)
    {
        if(!ptr) {
            return;
        }
        cudaFree(ptr);
        ptr = nullptr;
    }

    template <typename PtrT>
    inline
    PtrT* alloc_device_ptr(uint64_t size)
    {
        PtrT * ptr = nullptr;
        cudaMalloc(&ptr, sizeof(*ptr) * size);
        return ptr;
    }

    template <typename ValueT>
    inline
    void device_memset(ValueT * dst, ValueT value, std::size_t size)
    {
        cudaMemset(dst, value, sizeof(*dst) * size);
    }

    template <typename ValueT>
    inline
    void copy_from_device(ValueT * dst, ValueT const * src, std::size_t size)
    {
        cudaMemcpy(dst, src, sizeof(ValueT) * size, cudaMemcpyDeviceToHost);
    }

    template <typename ValueT>
    inline
    void copy_to_device(ValueT * dst, ValueT const* src, std::size_t size)
    {
        cudaMemcpy(dst, src, sizeof(ValueT) * size, cudaMemcpyHostToDevice);
    }

}