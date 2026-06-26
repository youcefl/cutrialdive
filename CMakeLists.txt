#
# MIT License
# Created on 2026.05.14
# Copyright (c) Youcef Lemsafer
# See LICENSE file for details.
#

cmake_minimum_required(VERSION 3.20)
project(cutrialdive VERSION 1.0.0 LANGUAGES CXX)
add_compile_definitions(CUTRIALDIVE_VERSION="${PROJECT_VERSION}")

option(CUTRIALDIVE_ENABLE_GPU "Build with CUDA GPU support" OFF)
option(CUTRIALDIVE_ENABLE_PRP "Enable PRP testing (requires GWNUM)" OFF)
option(CUTRIALDIVE_USE_LFP "Use LFP instead of primesieve" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_STANDARD_REQUIRED ON)


if(CUTRIALDIVE_ENABLE_GPU)
    set(CMAKE_CUDA_STANDARD 20)
    set(CUDA_STANDARD_REQUIRED ON)
    find_package(CUDAToolkit REQUIRED)
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES "75;80;86;89;90;100;103;120;121")
    endif()
    set(CUDA_ARCHS ${CMAKE_CUDA_ARCHITECTURES})
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

if(NOT GMP_INCLUDE_DIR OR NOT GMP_LIBRARY)
    message(FATAL_ERROR "GMP not found. Install libgmp-dev or set GMP_ROOT to the install prefix.")
endif()

if(CUTRIALDIVE_ENABLE_PRP)
    add_compile_definitions(CUTRIALDIVE_ENABLE_PRP)
endif()

if(CUTRIALDIVE_USE_LFP)
    add_compile_definitions(CUTRIALDIVE_HAS_LFP)
else()
    find_package(primesieve REQUIRED static)
    add_compile_definitions(CUTRIALDIVE_HAS_PRIMESIEVE)
endif()

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
list(FILTER SOURCES EXCLUDE REGEX "main\\.cpp$")

add_library(cutrialdive_core STATIC ${SOURCES})

add_executable(cutrialdive src/main.cpp)

if(CUTRIALDIVE_ENABLE_GPU)
    string(REPLACE ";" ", " CUDA_ARCHS_STR "${CMAKE_CUDA_ARCHITECTURES}")
    target_compile_definitions(cutrialdive_core PRIVATE CUDA_ARCHS="${CUDA_ARCHS_STR}")
    target_link_libraries(cutrialdive_core PRIVATE CUDA::cudart)
endif()

target_include_directories(cutrialdive_core PUBLIC ${GMP_INCLUDE_DIR})
target_include_directories(cutrialdive_core PUBLIC src)
if(CUTRIALDIVE_USE_LFP)
    target_include_directories(cutrialdive_core PUBLIC ${LFP_INCLUDE_DIR})
else()
    target_link_libraries(cutrialdive_core PUBLIC primesieve::primesieve)
endif()
target_link_libraries(cutrialdive PRIVATE cutrialdive_core)
target_link_libraries(cutrialdive PRIVATE ${GMP_LIBRARY})
if(CUTRIALDIVE_ENABLE_PRP)
    target_include_directories(cutrialdive_core PUBLIC ${GWNUM_INCLUDE})
    target_link_libraries(cutrialdive PRIVATE ${GWNUM_LIB})
endif()
if(CUTRIALDIVE_ENABLE_GPU)
    target_link_libraries(cutrialdive PRIVATE CUDA::cudart)
endif()

install(TARGETS cutrialdive RUNTIME DESTINATION bin)

enable_testing()

add_subdirectory(tests)
