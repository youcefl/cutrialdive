#
# MIT License
# Created on 2026.05.23
# Copyright (c) Youcef Lemsafer
# See LICENSE file for details
#

include(FetchContent)

FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.8.1
)

FetchContent_MakeAvailable(Catch2)

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

set(SOURCES
    checkpointing_tests.cpp
    filesystem_helper.cpp
    mersenne_tests.cpp
    proth_tests.cpp
    smarandache_tests.cpp
    smarandache_gp_tests.cpp
    test_modular.cpp
    tf_tests.cpp
)
if(CUTRIALDIVE_ENABLE_GPU)
    list(APPEND SOURCES ./tf.cu)
else()
    list(APPEND SOURCES ./tf.cpp)
endif()


add_executable(cutrialdive_tests ${SOURCES})

target_link_libraries(cutrialdive_tests PRIVATE cutrialdive_core)
target_link_libraries(cutrialdive_tests PRIVATE ${GMP_LIBRARY})
target_link_libraries(cutrialdive_tests PRIVATE
    Catch2::Catch2WithMain
)

include(CTest)
include(Catch)

catch_discover_tests(cutrialdive_tests)
