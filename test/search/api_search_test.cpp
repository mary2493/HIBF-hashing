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
    config.reads = data("query.fq");
    config.index_file = data("kmer.index");
    config.kmer_size = 20;
    config.window_size = 20;
    config.hash = hash_type::minimiser;

    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();

    EXPECT_NO_THROW(search(config));

    std::string const std_cout = testing::internal::GetCapturedStdout();
    std::string const std_cerr = testing::internal::GetCapturedStderr();

    std::string const expected_cout{"The following hits were found:\n"
                                    "query1: [0]\n"
                                    "query2: [1]\n"
                                    "query3: [2]\n"};

    EXPECT_EQ(expected_cout, std_cout);
    EXPECT_EQ("", std_cerr);
}

TEST_F(api_search_test, default_config_minimiser)
{
    configuration config{};
    config.reads = data("query.fq");
    config.index_file = data("minimiser.index");
    config.kmer_size = 20;
    config.window_size = 24;
    config.hash = hash_type::minimiser;

    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();

    EXPECT_NO_THROW(search(config));

    std::string const std_cout = testing::internal::GetCapturedStdout();
    std::string const std_cerr = testing::internal::GetCapturedStderr();

    std::string const expected_cout{"The following hits were found:\n"
                                    "query1: [0]\n"
                                    "query2: [1]\n"
                                    "query3: [2]\n"};

    EXPECT_EQ(expected_cout, std_cout);
    EXPECT_EQ("", std_cerr);
}

TEST_F(api_search_test, default_config_syncmer)
{
    configuration config{};
    config.reads = data("query.fq");
    config.index_file = data("syncmer.index");
    config.kmer_size = 15;
    config.s = 11;
    config.t = 2;
    config.hash = hash_type::syncmer;

    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();

    EXPECT_NO_THROW(search(config));

    std::string const std_cout = testing::internal::GetCapturedStdout();
    std::string const std_cerr = testing::internal::GetCapturedStderr();

    std::string const expected_cout{"The following hits were found:\n"
                                    "query1: [0]\n"
                                    "query2: [1]\n"
                                    "query3: [2]\n"};

    EXPECT_EQ(expected_cout, std_cout);
    EXPECT_EQ("", std_cerr);
}
