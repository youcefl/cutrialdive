#
# Created on 2026.05.14
# Copyright (c) Youcef Lemsafer
#

cmake_minimum_required(VERSION 3.20)
project(cutrialdive VERSION 0.0.1 LANGUAGES CXX)
add_compile_definitions(CUTRIALDIVE_VERSION="${PROJECT_VERSION}")

option(CUTRIALDIVE_ENABLE_GPU "Build with CUDA GPU support" OFF)
option(CUTRIALDIVE_ENABLE_PRP "Enable PRP testing (requires GWNUM)" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_STANDARD_REQUIRED ON)


if(CUTRIALDIVE_ENABLE_GPU)
    set(CMAKE_CUDA_STANDARD 20)
    set(CUDA_STANDARD_REQUIRED ON)
    find_package(CUDAToolkit REQUIRED)
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES "61;75;86;89;100;120-virtual")
    endif()
    enable_language(CUDA)
    add_compile_definitions(CUTRIALDIVE_ENABLE_GPU)
    add_compile_options($<$<COMPILE_LANGUAGE:CUDA>:-std=c++20>)
endif()

if(CUTRIALDIVE_ENABLE_PRP)
    enable_language(C)

    # GWNUM is required for PRP
    if(NOT DEFINED GWNUM_INCLUDE OR NOT DEFINED GWNUM_LIB)
        message(FATAL_ERROR "PRP enabled but GWNUM_INCLUDE and/or GWNUM_LIB not set.\n"
                            "Pass -DGWNUM_INCLUDE=/path/to/gwnum/include -DGWNUM_LIB=/path/to/gwnum.a")    
    endif()
endif()

find_path(GMP_INCLUDE_DIR gmp.h
    HINTS ${GMP_ROOT} $ENV{GMP_ROOT}
    PATH_SUFFIXES include
)
find_library(GMP_LIBRARY gmp
    HINTS ${GMP_ROOT} $ENV{GMP_ROOT}
    PATH_SUFFIXES lib
)
find_package(OpenMP REQUIRED)

if(NOT GMP_INCLUDE_DIR OR NOT GMP_LIBRARY)
    message(FATAL_ERROR "GMP not found. Install libgmp-dev or set GMP_ROOT to the install prefix.")
endif()

if(CUTRIALDIVE_ENABLE_PRP)
    add_compile_definitions(CUTRIALDIVE_ENABLE_PRP)
endif()

find_package(primesieve REQUIRED static)
add_compile_definitions(CUTRIALDIVE_HAS_PRIMESIEVE)

file(GLOB_RECURSE SOURCES src/*.cpp)
if(CUTRIALDIVE_ENABLE_GPU)
    file(GLOB_RECURSE CUDA_SOURCES src/*.cu)
    list(APPEND SOURCES ${CUDA_SOURCES})
endif()

if(CUTRIALDIVE_ENABLE_PRP)
    set_source_files_properties(src/gw_utility.c PROPERTIES LANGUAGE C)
    list(APPEND SOURCES src/gw_utility.c)
else()
    list(FILTER SOURCES EXCLUDE REGEX "prp\\.cpp$")
endif()

if(CUTRIALDIVE_ENABLE_GPU)
    list(FILTER SOURCES EXCLUDE REGEX "autotests\\.cpp$")
    list(FILTER SOURCES EXCLUDE REGEX "trial_factoring\\.cpp$")
endif()

add_executable(cutrialdive ${SOURCES})

target_include_directories(cutrialdive PRIVATE ${GMP_INCLUDE_DIR})
target_link_libraries(cutrialdive PRIVATE ${GMP_LIBRARY})
target_link_libraries(cutrialdive PRIVATE OpenMP::OpenMP_CXX)
target_include_directories(cutrialdive PRIVATE src)
target_link_libraries(cutrialdive PRIVATE primesieve::primesieve)
if(CUTRIALDIVE_ENABLE_PRP)
    target_include_directories(cutrialdive PRIVATE ${GWNUM_INCLUDE})
    target_link_libraries(cutrialdive PRIVATE ${GWNUM_LIB})
endif()
if(CUTRIALDIVE_ENABLE_GPU)
    target_link_libraries(cutrialdive PRIVATE CUDA::cudart)
endif()

install(TARGETS cutrialdive RUNTIME DESTINATION bin)


enable_testing()

add_subdirectory(tests)
