#include "search/search.hpp"

#include <iostream>
#include <ranges>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/search.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

#include <cereal/archives/binary.hpp>
#include <hibf/hierarchical_interleaved_bloom_filter.hpp>
#include <hibf/config.hpp>  
using namespace seqan3::literals;

void search(configuration const &)
{

    seqan::hibf::hierarchical_interleaved_bloom_filter hibf;
    try
    {
        std::ifstream is{config.index_file, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(hibf);
    }
    catch (std::exception const & e)
    {
        std::cerr << "Error: Failed to load index from " << config.index_file << ": " << e.what() << '\n';
        return;
    }

    seqan3::sequence_file_input<seqan3::sequence_file_input_default_traits_dna> reads_file{config.reads};
    std::vector<std::string> matched_reads;

    seqan3::debug_stream << "=====   Running on a single text   =====\n"
                         << "The following hits were found:\n";

    auto agent = hibf.membership_agent();
    size_t const threshold = 1u;

    /*  auto & result1 = agent.membership_for(query1, 2u);

    // query1 hits in user_bin_1 and user_bin_3, which have the IDs 0 and 2, respectively.
    for (uint64_t hit_user_bin : result1)
        std::cout << hit_user_bin << ' '; // The results are not sorted: 2 0
    std::cout << '\n';
    
    for (auto & record : reads_fin)
{
    auto & result = agent.membership_for(record.sequence() | minimiser_view, threshold);
    seqan3::debug_stream << read.id() << ": " << result << '\n';
}
    */

    for (auto & record : reads_file)
    {
        auto window_size = config.kmer_size;
        auto kmer_view = record.sequence() | seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{config.kmer_size}}, seqan3::window_size{window_size});
        auto & result = agent.membership_for(kmer_view, threshold);
        seqan3::debug_stream << record.id() << ": " << result << '\n';
    }
}