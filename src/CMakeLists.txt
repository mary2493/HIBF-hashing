# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

cmake_minimum_required (VERSION 3.25)

# An object library (without main) to be used in multiple targets.
# You can add more external include paths of other projects that are needed for your project.
add_library (HIBF-hashing_lib STATIC build/build.cpp build/run_build.cpp search/search.cpp search/run_search.cpp)
target_include_directories (HIBF-hashing_lib PUBLIC "${HIBF-hashing_SOURCE_DIR}/include")
target_link_libraries (HIBF-hashing_lib PUBLIC seqan3::seqan3 sharg::sharg seqan::hibf)

target_compile_options (HIBF-hashing_lib PUBLIC "-pedantic" "-Wall" "-Wextra")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # Disable warning about std::hardware_destructive_interference_size not being ABI-stable.
    if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 12)
        target_compile_options (HIBF-hashing_lib PUBLIC "-Wno-interference-size")
    endif ()

    # Warn about failed return value optimization.
    if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 14)
        target_compile_options (HIBF-hashing_lib PUBLIC "-Wnrvo")
    endif ()
endif ()

# Add the application.
add_executable (HIBF-hashing main.cpp)
target_link_libraries (HIBF-hashing PRIVATE HIBF-hashing_lib)
