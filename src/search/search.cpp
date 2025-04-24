// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include "search/search.hpp"

#include <iostream>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/search.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

#include "contrib/syncmer.hpp"
#include "dna4_traits.hpp"
#include "index_data.hpp"
#include <cereal/archives/binary.hpp>
#include <hibf/config.hpp>
#include <hibf/hierarchical_interleaved_bloom_filter.hpp>
#include <threshold/threshold.hpp>

void search(configuration const & config)
{
    using sequence_file_t = seqan3::sequence_file_input<dna4_traits>;

    myindex index{};
    index.load(config.index_file);

    if (index.hash == hash_type::syncmer)
    {
        if (config.error > 0)
            throw std::runtime_error{"Syncmer does not support errors."};
    }

    threshold::threshold const thresholder = [&]() -> threshold::threshold
    {
        size_t const first_sequence_size = [&]()
        {
            sequence_file_t fin{config.reads};
            auto & record = *fin.begin();
            return record.sequence().size();
        }();

        return {threshold::threshold_parameters{.window_size = index.window_size,
                                                .shape = seqan3::ungapped{index.kmer_size},
                                                .query_length = first_sequence_size,
                                                .errors = config.error}};
    }();

    sequence_file_t reads_file{config.reads};

    auto agent = index.hibf.membership_agent();

    //for storing the results
    std::vector<std::string> results;
    std::string current_read{};
    std::vector<uint64_t> hashes;

    auto get_results = [&](auto & record)
    {
        auto & result = agent.membership_for(hashes, thresholder.get(hashes.size()));
        agent.sort_results();

        current_read.clear();
        current_read += record.id() + ": [";

        for (size_t i = 0; i < result.size(); ++i)
        {
            current_read += std::to_string(result[i]);
            if (i < result.size() - 1)
                current_read += ",";
        }
        current_read += "]\n";

        // store the result in the vector
        results.push_back(current_read);
    };

    if (index.hash != hash_type::syncmer)
    {
        auto hash_adaptor =
            seqan3::views::minimiser_hash(seqan3::ungapped{index.kmer_size}, seqan3::window_size{index.window_size});

        for (auto & record : reads_file)
        {
            if (record.sequence().size() < index.window_size)
            {
                throw std::runtime_error{"read in file is shorter than the k-mer/window size."};
            }

            auto minimiser_view = record.sequence() | hash_adaptor | std::views::common;
            hashes.clear();
            hashes.assign(minimiser_view.begin(), minimiser_view.end());
            get_results(record);
        }
    }

    else
    {
        auto hash_adaptor =
            seqan3::views::syncmer({.kmer_size = index.kmer_size, .smer_size = index.s, .offset = index.t});

        for (auto & record : reads_file)
        {
            if (record.sequence().size() < index.kmer_size)
            {
                throw std::runtime_error{"read in file is shorter than the k-mer size."};
            }
            auto syncmer_view = record.sequence() | hash_adaptor | std::views::common;
            hashes.clear();
            hashes.assign(syncmer_view.begin(), syncmer_view.end());
            get_results(record);
        }
    }

    std::cout << "The following hits were found:\n";
    //print to console
    for (auto & record : results)
    {
        std::cout << record;
    }
    //save to file
    std::ofstream result_out{config.search_output};
    for (auto & record : results)
    {
        result_out << record;
    }
}
