// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <gtest/gtest.h>

#include "../app_test.hpp"
#include <search/search.hpp>

// To prevent issues when running multiple API tests in parallel, give each API test unique names:
struct api_search_test : public app_test
{};

// Test für den `search`-Prozess mit einer Standardkonfiguration:
TEST_F(api_search_test, default_config)
{
    configuration config{};
    config.reads = data("test_reads.fasta");
    config.index_file = data("test_index.bin");
    config.kmer_size = 20;
    config.threshold = 1; 

    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();

    search(config);

    std::string const std_cout = testing::internal::GetCapturedStdout();
    std::string const std_cerr = testing::internal::GetCapturedStderr();

    // Erwartete Standardausgabe und Fehlerausgabe
    std::string const expected_cout{"==============================\n"
                                    "The following hits were found:\n"
                                    "read1: [0,1]\n"
                                    "read2: []\n"
                                    "read3: [0]\n"
                                    "read4: []\n"};
    std::string const expected_cerr{}; 

    EXPECT_EQ(expected_cout, std_cout);
    EXPECT_EQ(expected_cerr, std_cerr);
}