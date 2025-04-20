// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include "build/run_build.hpp"

#include "build/build.hpp"
#include "configuration.hpp"

// This tells sharg how to map strings to enum values and vice versa.
auto enumeration_names(hash_type)
{
    return std::unordered_map<std::string_view, hash_type>{{"kmer", hash_type::kmer},
                                                           {"minimiser", hash_type::minimiser},
                                                           {"syncmer", hash_type::syncmer}};
}

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
                      .validator = sharg::output_file_validator{sharg::output_file_open_options::open_or_create}});

    parser.add_option(config.kmer_size,
                      sharg::config{.short_id = 'k',
                                    .long_id = "kmer",
                                    .description = "The k-mer size to use.",
                                    .validator = sharg::arithmetic_range_validator{1, 32}});

    parser.add_option(config.window_size,
                      sharg::config{.short_id = 'w',
                                    .long_id = "window",
                                    .description = "The window size for minimisers.",
                                    .default_message = "k-mer size",
                                    .validator = sharg::arithmetic_range_validator{1, 200}});

    parser.add_option(config.window_size,
                      sharg::config{.short_id = 's',
                                    .description = "length of the smaller s-mer",
                                    .validator = sharg::arithmetic_range_validator{1, 32}});

    parser.add_option(config.window_size,
                      sharg::config{.short_id = 't',
                                    .description = "position within the k-mer at which the minimal s-mer must occur "
                                                   "for the k-mer to be selected as a syncmer",
                                    .validator = sharg::arithmetic_range_validator{1, 32}});

    parser.add_option(config.hash,
                      sharg::config{.short_id = 'm',
                                    .long_id = "mode",
                                    .description = "Hash type to use.",
                                    .validator = sharg::value_list_validator{
                                        (sharg::enumeration_names<hash_type> | std::views::values)}});

    try
    {
        parser.parse(); // Trigger command line parsing.

        // Make window default to kmer size if not set.
        if (!parser.is_option_set("window"))
        {
            config.window_size = config.kmer_size;
        }
        else if (config.window_size < config.kmer_size)
        {
            throw std::runtime_error{"Window size must be greater than or equal to k-mer size."};
        }
        else // If window provided, user wants minimiser.
        {
            config.hash = hash_type::minimiser;
        }
    }
    catch (sharg::parser_error const & ext) // Catch user errors.
    {
        std::cerr << "Parsing error. " << ext.what() << '\n'; // Give error message.
        std::exit(-1);
    }

    build(config);
}
