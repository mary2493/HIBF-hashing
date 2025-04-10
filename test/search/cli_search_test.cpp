// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include "../app_test.hpp"

// To prevent issues when running multiple CLI tests in parallel, give each CLI test unique names:
struct cli_search_test : public app_test
{};

TEST_F(cli_search_test, no_options)
{
    app_test_result const result = execute_app("HIBF-hashing", "search");
    std::string_view const expected{"HIBF-hashing-search\n"
                                    "===================\n"
                                    "    Try -h or --help for more information.\n"};

    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, "");
}

TEST_F(cli_search_test, missing_required_argument)
{
    app_test_result const result = execute_app("HIBF-hashing", "search", "-e", "-o", "w", "-t");
    std::string_view const expected{"Parsing error. Option -i/--index is required but not set.\n"};

    EXPECT_FAILURE(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, expected);
}

TEST_F(cli_search_test, with_arguments)
{
    app_test_result const result = execute_app("HIBF-hashing",
                                               "search",
                                               "--index",
                                               data("test_index_kmer.bin"),
                                               "--reads",
                                               data("reads.fasta"),
                                               "--error 1",
                                               "--output result.out",
                                               "--window 20",
                                               "--type kmer");

    std::string const expected{"The following hits were found:\n"
                               "read1: [1,2]\n"
                               "read2: []\n"
                               "read3: [1]\n"

    };

    EXPECT_SUCCESS(result);
    EXPECT_EQ(result.out, expected);
}

TEST_F(cli_search_test, missing_path)
{
    app_test_result const result = execute_app("HIBF-hashing",
                                               "search",
                                               "--index",
                                               data("test_index_kmer.bin"),
                                               "--reads",
                                               data("reads.fasta"),
                                               "--error 1",
                                               "--output",
                                               data("search_test.txt"),
                                               "-o",
                                               "");

    EXPECT_FAILURE(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err, "Parsing error. Missing value for option -o\n");
}

TEST_F(cli_search_test, invalid_input_path)
{
    app_test_result const result = execute_app("HIBF-hashing", "search", "--index does/not/exist");

    EXPECT_FAILURE(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err,
              "Parsing error. Validation failed for option -i/--index: The file \"does/not/exist\" does not exist!\n");
}

TEST_F(cli_search_test, invalid_output_path)
{
    app_test_result const result = execute_app("HIBF-hashing",
                                               "search",
                                               "--index",
                                               data("test_index_kmer.bin"),
                                               "--reads",
                                               data("reads.fasta"),
                                               "--error 1",
                                               "--output does/not/exist");

    EXPECT_FAILURE(result);
    EXPECT_EQ(result.out, "");
    EXPECT_EQ(result.err,
              "Parsing error. Validation failed for option -o/--output: Cannot write \"does/not/exist\"!\n");
}
