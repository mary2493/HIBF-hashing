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
#include <hibf/config.hpp>
#include <hibf/hierarchical_interleaved_bloom_filter.hpp>
using namespace seqan3::literals;

void search(configuration const & config)
{

    seqan::hibf::hierarchical_interleaved_bloom_filter hibf;

    std::ifstream is{config.index_file, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(hibf);

    seqan3::sequence_file_input<seqan3::sequence_file_input_default_traits_dna> reads_file{config.reads};
    std::vector<std::string> matched_reads;

    auto agent = hibf.membership_agent();
    size_t threshold = config.threshold;

    std::ofstream result_out{config.search_output};
    std::cout << "The following hits were found:\n";
    for (auto & record : reads_file)
    {
        auto kmer_view =
            record.sequence() | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{config.kmer_size}});
        auto & result = agent.membership_for(kmer_view, threshold);

        //output console
        std::cout << record.id() << ": [";
        for (size_t i = 0; i < result.size(); ++i)
        {
            std::cout << result[i];
            if (i < result.size() - 1)
                std::cout << ",";
        }
        std::cout << "]\n";

        //save to file
        result_out << record.id() << ": ";
        for (size_t i = 0; i < result.size(); ++i)
        {
            result_out << result[i];
            if (i < result.size() - 1)
                result_out << ", ";
        }
        result_out << '\n';
    }
}
