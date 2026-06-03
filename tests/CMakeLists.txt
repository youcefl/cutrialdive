#
# Created on 2026.05.23
# Copyright (c) Youcef Lemsafer
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

add_executable(cutrialdive_tests
    test_modular.cpp
    mersenne_tests.cpp
    smarandache_tests.cpp
    smarandache_gp_tests.cpp
    ../src/number_sequences/smarandache.cpp
)

target_include_directories(cutrialdive_tests PRIVATE ${GMP_INCLUDE_DIR})
target_include_directories(cutrialdive_tests PRIVATE ${PROJECT_SOURCE_DIR}/src)

target_link_libraries(cutrialdive_tests PRIVATE ${GMP_LIBRARY})
target_link_libraries(cutrialdive_tests PRIVATE
    Catch2::Catch2WithMain
)

include(CTest)
include(Catch)

catch_discover_tests(cutrialdive_tests)
