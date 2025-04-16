// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <gtest/gtest.h>

#include "../app_test.hpp"
#include <build/build.hpp>

// To prevent issues when running multiple API tests in parallel, give each API test unique names:
struct api_build_test : public app_test
{};

TEST_F(api_build_test, default_config_kmer)
{
    configuration config{};
    config.file_list_path = data("file_list_for_tests.txt");
    config.index_output = "new_kmer.index";
    config.kmer_size = 20;
    config.hash = hash_type::kmer;
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();

    EXPECT_NO_THROW(build(config));

    std::string const std_cout = testing::internal::GetCapturedStdout();
    std::string const std_cerr = testing::internal::GetCapturedStderr();

    std::string const expected_cout{"HIBF index built and saved to \"new_kmer.index\"\n"
                                    "Successfully processed 3 files.\n"};

    EXPECT_EQ(expected_cout, std_cout);
    EXPECT_EQ("", std_cerr);
}

TEST_F(api_build_test, default_config_minimiser)
{
    configuration config{};
    config.file_list_path = data("file_list_for_tests.txt");
    config.index_output = "new_minimiser.index";
    config.kmer_size = 18;
    config.window_size = 20;
    config.hash = hash_type::minimiser;

    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();

    EXPECT_NO_THROW(build(config));

    std::string const std_cout = testing::internal::GetCapturedStdout();
    std::string const std_cerr = testing::internal::GetCapturedStderr();

    std::string const expected_cout{"HIBF index built and saved to \"new_minimiser.index\"\n"
                                    "Successfully processed 3 files.\n"};

    EXPECT_EQ(expected_cout, std_cout);
    EXPECT_EQ("", std_cerr);
}
