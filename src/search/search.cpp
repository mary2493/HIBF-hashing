// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include "search/search.hpp"

#include <iostream>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

#include "contrib/syncmer.hpp"
#include "dna4_traits.hpp"
#include "index_data.hpp"
#include <cereal/archives/binary.hpp>
#include <hibf/config.hpp>
#include <hibf/hierarchical_interleaved_bloom_filter.hpp>
#include <threshold/threshold.hpp>

threshold::threshold get_thresholder(configuration const & config, myindex const & index)
{
    size_t const first_sequence_size = [&]()
    {
        seqan3::sequence_file_input<dna4_traits> fin{config.reads};
        auto & record = *fin.begin();
        return record.sequence().size();
    }();

    return {threshold::threshold_parameters{.window_size = index.window_size,
                                            .shape = seqan3::ungapped{index.kmer_size},
                                            .query_length = first_sequence_size,
                                            .errors = config.error}};
}

void search(configuration const & config)
{
    myindex index{};
    index.load(config.index_file);
    auto agent = index.hibf.membership_agent();

    std::vector<std::string> results;
    std::string result_line{};
    std::array<char, std::numeric_limits<uint64_t>::digits10 + 1> buffer{};
    std::vector<uint64_t> hashes;

    threshold::threshold const thresholder = get_thresholder(config, index);
    auto get_results = [&](auto & record)
    {
        auto & result = agent.membership_for(hashes, thresholder.get(hashes.size()));
        agent.sort_results();

        result_line.clear();
        result_line += record.id() + ": [";

        for (auto && bin : result)
        {
            auto conv = std::to_chars(buffer.data(), buffer.data() + buffer.size(), bin);
            assert(conv.ec == std::errc{});
            std::string_view sv{buffer.data(), conv.ptr};
            result_line += sv;
            result_line += ',';
        }

        if (result_line.back() == ',')
            result_line.pop_back();

        result_line += "]\n";

        // store the result in the vector
        results.push_back(result_line);
    };

    seqan3::sequence_file_input<dna4_traits> reads_file{config.reads};
    auto process = [&](auto hash_adaptor)
    {
        for (auto & record : reads_file)
        {
            auto view = record.sequence() | hash_adaptor | std::views::common;
            hashes.clear();
            hashes.assign(view.begin(), view.end());
            get_results(record);
        }
    };

    if (index.hash != hash_type::syncmer)
    {
        process(
            seqan3::views::minimiser_hash(seqan3::ungapped{index.kmer_size}, seqan3::window_size{index.window_size}));
    }
    else
    {
        process(seqan3::views::syncmer({.kmer_size = index.kmer_size, .smer_size = index.s, .offset = index.t}));
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
