// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include "build/build.hpp"

#include <algorithm> // for std::all_of
#include <cctype>    // for isspace
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <sharg/validators.hpp>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

#include "contrib/syncmer.hpp"
#include "dna4_traits.hpp"
#include "index_data.hpp"
#include <cereal/archives/binary.hpp>
#include <hibf/config.hpp>
#include <hibf/hierarchical_interleaved_bloom_filter.hpp>

template <hash_type hash>
std::function<void(size_t, seqan::hibf::insert_iterator &&)>
get_input_fn_impl(configuration const & config, std::vector<std::string> const & user_bin_paths)
{
    using sequence_file_t = seqan3::sequence_file_input<dna4_traits>;

    auto view = [&]()
    {
        if constexpr (hash == hash_type::minimiser)
        {
            return seqan3::views::minimiser_hash(seqan3::ungapped{config.kmer_size},
                                                 seqan3::window_size{config.window_size});
        }
        else
        {
            static_assert(hash == hash_type::syncmer);
            return seqan3::views::syncmer({.kmer_size = config.kmer_size, .smer_size = config.s, .offset = config.t});
        }
    }();

    return [&, view](size_t const user_bin_id, seqan::hibf::insert_iterator it)
    {
        sequence_file_t fin{user_bin_paths[user_bin_id]};
        for (auto & record : fin)
        {
            std::ranges::copy(record.sequence() | view, it);
        }
    };
}

std::function<void(size_t, seqan::hibf::insert_iterator &&)>
get_input_fn(configuration const & config, std::vector<std::string> const & user_bin_paths)
{
    switch (config.hash)
    {
    case hash_type::minimiser:
        return get_input_fn_impl<hash_type::minimiser>(config, user_bin_paths);
    case hash_type::syncmer:
        return get_input_fn_impl<hash_type::syncmer>(config, user_bin_paths);
    default:
        throw std::runtime_error{"Invalid hash type."};
    }
}

std::vector<std::string> parse_user_bins(std::filesystem::path const & file)
{
    std::vector<std::string> user_bin_paths;
    std::ifstream file_list{file};
    std::string current_line;
    sharg::input_file_validator fasta_validator{{"fasta", "fa", "fna"}};

    //Each FASTA file is opened, and the k-mers are extracted from it.
    //These kmers are stored in all_bins_together, with each file corresponding to a "User Bin" in the HIBF
    while (std::getline(file_list, current_line))
    {

        /* Checking whether the line is empty or consists only of whitespace
        .empty() checks if the string is empty, and std::all_of checks if all characters in the string are whitespace. It returns true
        if every character is a whitespace (e.g., space, tab, newline) and false if at least one character is not a whitespace */
        if (current_line.empty()
            || std::ranges::all_of(current_line,
                                   [](unsigned char c)
                                   {
                                       return std::isspace(c);
                                   }))
        {
            throw std::runtime_error{"Empty line or invalid entry in the file list."};
        }

        fasta_validator(current_line);
        user_bin_paths.push_back(current_line);
    }

    if (user_bin_paths.empty())
        throw std::runtime_error{"No valid files found in the file list."};

    return user_bin_paths;
}

void build(configuration const & config)
{
    std::vector<std::string> const user_bin_paths = parse_user_bins(config.file_list_path);

    auto input_fn = get_input_fn(config, user_bin_paths);

    seqan::hibf::config hibf_config{.input_fn = input_fn,                         // required
                                    .number_of_user_bins = user_bin_paths.size(), // required
                                    .number_of_hash_functions = 2u,
                                    .maximum_fpr = 0.05,
                                    .threads = 1u};

    // The HIBF constructor will determine a hierarchical layout for the user bins and build the filter
    seqan::hibf::hierarchical_interleaved_bloom_filter hibf{hibf_config};

    //The indices can also be stored and loaded from disk by using cereal
    myindex index{config.kmer_size, std::move(hibf), config.window_size, config.s, config.t, config.hash};
    index.store(config.index_output);

    std::cout << "HIBF index built and saved to " << config.index_output << "\n";
    std::cout << "Successfully processed " << user_bin_paths.size() << " files.\n";
}
