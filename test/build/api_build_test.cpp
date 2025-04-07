// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <gtest/gtest.h>

#include "../app_test.hpp"
#include <build/build.hpp>

// To prevent issues when running multiple API tests in parallel, give each API test unique names:
struct api_build_test : public app_test
{};

TEST_F(api_build_test, default_config)
{
    configuration config{};
    config.file_list_path = data("file_list_for_tests.txt");
    config.index_output = "new.index";
    config.kmer_size = 20;
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();

    build(config);

    std::string const std_cout = testing::internal::GetCapturedStdout();
    std::string const std_cerr = testing::internal::GetCapturedStderr();

    std::string const expected_cout{"HIBF index built and saved to \"new.index\"\n"
                                    "Successfully processed 2 files.\n"};
    std::string const expected_cerr{"Error: Could not parse file " + data("file_test1.fasta").string() + ".\n"
                                    "Error: Sequence in file " + data("file_test2.fasta").string() + " is shorter than the k-mer size. Skipping sequence.\n"
                                    "Error: Empty line or invalid entry in the file list.\n"
                                    "Error: Unsupported file format for file " + data("file_test4.txt").string() + ".\n"};

    EXPECT_EQ(expected_cout, std_cout);
    EXPECT_EQ(expected_cerr, std_cerr);
}
