/*
* Creation date: 2026.05.04
* Copyright (c) Youcef Lemsafer
*/
#pragma once

#if defined __CUDA_ARCH__
#  define CUTRIALDIVE_IS_CUDA 1
#  define CUTRIALDIVE_DEVICE __device__
#  define CUTRIALDIVE_HOST __host__
#  define CUTRIALDIVE_DEVICE_AND_HOST CUTRIALDIVE_DEVICE CUTRIALDIVE_HOST
#else
#  define CUTRIALDIVE_IS_CUDA 0
#  define CUTRIALDIVE_DEVICE
#  define CUTRIALDIVE_HOST
#  define CUTRIALDIVE_DEVICE_AND_HOST
#endif
