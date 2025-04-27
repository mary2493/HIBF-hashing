// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <sharg/all.hpp>

#include "build/run_build.hpp"
#include "search/run_search.hpp"

int main(int argc, char ** argv)
{
    // Parser
    sharg::parser parser{"HIBF-hashing",
                         argc,
                         argv,
                         sharg::update_notifications{sharg::update_notifications::off},
                         {"build", "search"}};

    // General information.
    parser.info.author = "Mariya";
    parser.info.version = "1.0.0";

    try
    {
        parser.parse();

        // hold a reference to the sub_parser
        sharg::parser & sub_parser = parser.get_sub_parser();

        if (sub_parser.info.app_name == std::string_view{"HIBF-hashing-build"})
            run_build(sub_parser);
        else if (sub_parser.info.app_name == std::string_view{"HIBF-hashing-search"})
            run_search(sub_parser);
    }
    catch (std::exception const & ext)
    {
        std::cerr << "[Error] " << ext.what() << '\n';
        return -1;
    }
    return 0;
}
