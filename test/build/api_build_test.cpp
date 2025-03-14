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
    testing::internal::CaptureStdout();
    build(config);
    std::string const std_cout = testing::internal::GetCapturedStdout();

    std::string const expected{"\"\"\n"
                               "\"index\"\n"
                               "20\n"};

    EXPECT_EQ(expected, std_cout);
}
