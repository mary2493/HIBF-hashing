// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <gtest/gtest.h>

#include "../app_test.hpp"
#include <search/search.hpp>

// To prevent issues when running multiple API tests in parallel, give each API test unique names:
struct api_search_test : public app_test
{};

TEST_F(api_search_test, default_config_kmer)
{
    configuration config{};
    config.reads = data("reads.fasta");
    config.index_file = data("test_index_kmer.bin");
    config.kmer_size = 20;
    config.hash = hash_type::kmer;

    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();

    EXPECT_NO_THROW(search(config));

    std::string const std_cout = testing::internal::GetCapturedStdout();
    std::string const std_cerr = testing::internal::GetCapturedStderr();

    std::string const expected_cout{"The following hits were found:\n"
                                    "read1: [1,2]\n"
                                    "read2: []\n"
                                    "read3: [1]\n"};

    EXPECT_EQ(expected_cout, std_cout);
    EXPECT_EQ("", std_cerr);
}

TEST_F(api_search_test, default_config_minimiser)
{
    configuration config{};
    config.reads = data("reads.fasta");
    config.index_file = data("test_index_minimiser.bin");
    config.kmer_size = 18;
    config.window_size = 20;
    config.hash = hash_type::minimiser;

    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();

    EXPECT_NO_THROW(search(config));

    std::string const std_cout = testing::internal::GetCapturedStdout();
    std::string const std_cerr = testing::internal::GetCapturedStderr();

    std::string const expected_cout{"The following hits were found:\n"
                                    "read1: [1,2]\n"
                                    "read2: []\n"
                                    "read3: [1]\n"};

    EXPECT_EQ(expected_cout, std_cout);
    EXPECT_EQ("", std_cerr);
}
