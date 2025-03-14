// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include "build/run_build.hpp"

#include "build/build.hpp"
#include "configuration.hpp"

void run_build(sharg::parser & parser)
{
    configuration config{};

    parser.add_option(config.file_list_path,
                      sharg::config{.short_id = 'i',
                                    .long_id = "input",
                                    .description = "A file containing one sequence file per line",
                                    .required = true,
                                    .validator = sharg::input_file_validator{}});

    parser.add_option(
        config.index_output,
        sharg::config{.short_id = 'o',
                      .long_id = "output",
                      .description = "Where to store the index.",
                      .validator = sharg::output_file_validator{sharg::output_file_open_options::create_new}});

    parser.add_option(config.kmer_size,
                      sharg::config{.short_id = 'k',
                                    .long_id = "kmer",
                                    .description = "The k-mer size to use.",
                                    .validator = sharg::arithmetic_range_validator{1, 32}});

    try
    {
        parser.parse(); // Trigger command line parsing.
    }
    catch (sharg::parser_error const & ext) // Catch user errors.
    {
        std::cerr << "Parsing error. " << ext.what() << '\n'; // Give error message.
        std::exit(-1);
    }

    build(config);
}
