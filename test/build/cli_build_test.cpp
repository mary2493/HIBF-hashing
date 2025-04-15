// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include "../app_test.hpp"

// To prevent issues when running multiple CLI tests in parallel, give each CLI test unique names:
struct cli_build_test : public app_test
{};

TEST_F(cli_build_test, no_options)
{
    app_test_result const result = execute_app("HIBF-hashing", "build");
    std::string_view const expected{"HIBF-hashing-build\n"
                                    "==================\n"
                                    "    Try -h or --help for more information.\n"};

    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, "");
}

TEST_F(cli_build_test, missing_required_argument)
{
    app_test_result const result = execute_app("HIBF-hashing", "build", "--kmer 20");
    std::string_view const expected{"Parsing error. Option -i/--input is required but not set.\n"};

    EXPECT_FAILURE(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, expected);
}

TEST_F(cli_build_test, with_arguments_kmer)
{
    app_test_result const result = execute_app("HIBF-hashing",
                                               "build",
                                               "--input",
                                               data("file_list_for_tests.txt"),
                                               "--output new_kmer.index",
                                               "--kmer 20",
                                               "--window 20",
                                               "--mode kmer");

    std::string const expected{"HIBF index built and saved to \"new_kmer.index\"\n"
                               "Successfully processed 3 files.\n"};

    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, expected);
}

TEST_F(cli_build_test, with_arguments_minimiser)
{
    app_test_result const result = execute_app("HIBF-hashing",
                                               "build",
                                               "--input",
                                               data("file_list_for_tests.txt"),
                                               "--output new_minimiser.index",
                                               "--kmer 20",
                                               "--window 22",
                                               "--mode minimiser");

    std::string const expected{"HIBF index built and saved to \"new_minimiser.index\"\n"
                               "Successfully processed 3 files.\n"};

    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, expected);
}

TEST_F(cli_build_test, missing_path)
{
    app_test_result const result =
        execute_app("HIBF-hashing", "build", "--input", data("file_list_for_tests.txt"), "-o", "");

    EXPECT_FAILURE(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, "Parsing error. Missing value for option -o\n");
}

TEST_F(cli_build_test, invalid_input_path)
{
    app_test_result const result = execute_app("HIBF-hashing", "build", "--input does/not/exist");

    EXPECT_FAILURE(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err,
              "Parsing error. Validation failed for option -i/--input: The file \"does/not/exist\" does not exist!\n");
}

TEST_F(cli_build_test, invalid_output_path)
{
    app_test_result const result =
        execute_app("HIBF-hashing", "build", "--input", data("file_list_for_tests.txt"), "--output does/not/exist");

    EXPECT_FAILURE(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err,
              "Parsing error. Validation failed for option -o/--output: Cannot write \"does/not/exist\"!\n");
}
