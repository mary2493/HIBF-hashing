// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include "build/run_build.hpp"

#include "build/build.hpp"
#include "configuration.hpp"

void add_shared_options(sharg::parser & parser, configuration & config)
{
    parser.add_subsection("General options");
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
}

void run_minimiser(sharg::parser & parser)
{
    configuration config{.hash = hash_type::minimiser};

    add_shared_options(parser, config);

    parser.add_subsection("Minimizer options");
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

    parser.parse();

    // Make window default to kmer size if not set.
    if (!parser.is_option_set("window"))
        config.window_size = config.kmer_size;
    else if (config.window_size < config.kmer_size)
        throw std::runtime_error{"Window size must be greater than or equal to k-mer size."};

    build(config);
}

void run_syncmer(sharg::parser & parser)
{
    configuration config{.hash = hash_type::syncmer};

    add_shared_options(parser, config);

    parser.add_subsection("Syncmer options");
    parser.add_option(config.kmer_size,
                      sharg::config{.short_id = 'k',
                                    .long_id = "kmer",
                                    .description = "The k-mer size to use.",
                                    .validator = sharg::arithmetic_range_validator{1, 32}});
    parser.add_option(config.s,
                      sharg::config{.short_id = 's',
                                    .long_id = "syncmer_s",
                                    .description = "length of the smaller s-mer.",
                                    .validator = sharg::arithmetic_range_validator{1, 32}});
    parser.add_option(config.t,
                      sharg::config{.short_id = 't',
                                    .long_id = "syncmer_t",
                                    .description = "position within the k-mer at which the minimal s-mer must occur.",
                                    .validator = sharg::arithmetic_range_validator{0, 32}});

    parser.parse();

    if (config.s >= config.kmer_size)
        throw std::invalid_argument{"Syncmer s-mer size must be smaller than k-mer size."};
    if (config.t > config.kmer_size - config.s)
        throw std::invalid_argument{"Syncmer offset t is out of bounds."};

    build(config);
}

void run_build(sharg::parser & parser)
{
    parser.add_subcommands({"minimiser", "syncmer"});
    parser.parse();

    sharg::parser & sub_parser = parser.get_sub_parser();

    if (sub_parser.info.app_name == std::string_view{"HIBF-hashing-build-minimiser"})
        run_minimiser(sub_parser);
    else if (sub_parser.info.app_name == std::string_view{"HIBF-hashing-build-syncmer"})
        run_syncmer(sub_parser);
}
