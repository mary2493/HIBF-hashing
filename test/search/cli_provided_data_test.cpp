// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include "../app_test.hpp"

// To prevent issues when running multiple CLI tests in parallel, give each CLI test unique names:
struct search_test : public app_test, public testing::WithParamInterface<std::tuple<std::string_view, uint16_t>>
{
    static std::filesystem::path get_index(std::string_view const hash_type)
    {
        if (hash_type == "kmer")
            return data("20_20.index");
        else if (hash_type == "minimiser")
            return data("20_24.index");
        else
            throw std::runtime_error{"Unknown hash type"};
    }

    static std::string_view get_expected_hits(std::string_view const hash_type, uint16_t const errors)
    {
        if (errors == 2u)
        {
            if (hash_type == "kmer")
                return "The following hits were found:\nquery1: [0,1,2,3]\nquery2: [0,1,2,3]\nquery3: [1,2,3]\n";
            else
                return "The following hits were found:\nquery1: [0]\nquery2: [1,3]\nquery3: [2]\n";
        }

        return "The following hits were found:\nquery1: [0]\nquery2: [1]\nquery3: [2]\n";
    }
};

TEST_P(search_test, run)
{
    auto const [hash_type, number_of_errors] = GetParam();

    app_test_result const result = execute_app("HIBF-hashing",
                                               "search",
                                               "--index",
                                               get_index(hash_type),
                                               "--error",
                                               number_of_errors,
                                               "--reads",
                                               data("query.fq"),
                                               "--output result.out");

    EXPECT_SUCCESS(result);
    ASSERT_TRUE(std::filesystem::exists("result.out"));

    EXPECT_EQ(result.out, get_expected_hits(hash_type, number_of_errors));
    EXPECT_EQ(result.err, "");
}

INSTANTIATE_TEST_SUITE_P(provided_data,
                         search_test,
                         testing::Combine(testing::Values("kmer", "minimiser"), testing::Values(0, 1, 2)),
                         [](testing::TestParamInfo<search_test::ParamType> const & info)
                         {
                             std::string name{std::get<0>(info.param)};
                             name += "_";
                             name += std::to_string(std::get<1>(info.param));
                             name += "_errors";
                             return name;
                         });
