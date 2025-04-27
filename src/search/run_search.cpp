// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include "search/run_search.hpp"

#include "configuration.hpp"
#include "search/search.hpp"

void run_search(sharg::parser & parser)
{
    configuration config{};

    parser.add_option(config.index_file,
                      sharg::config{.short_id = 'i',
                                    .long_id = "index",
                                    .description = "HIBF index file to load.",
                                    .required = true,
                                    .validator = sharg::input_file_validator{}});

    parser.add_option(config.reads,
                      sharg::config{.short_id = 'r',
                                    .long_id = "reads",
                                    .description = "reads to search for.",
                                    .required = true,
                                    .validator = sharg::input_file_validator{}});

    parser.add_option(config.error,
                      sharg::config{.short_id = 'e',
                                    .long_id = "error",
                                    .description = "The maximum number of errors allowed.",
                                    .validator = sharg::arithmetic_range_validator{0, 5}});

    parser.add_option(
        config.search_output,
        sharg::config{.short_id = 'o',
                      .long_id = "output",
                      .description = ".txt file to write the search results to.",
                      .validator = sharg::output_file_validator{sharg::output_file_open_options::open_or_create}});

    parser.parse();

    search(config);
}
