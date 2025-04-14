#include "search/search.hpp"

#include <iostream>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/search.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

#include "hash_utils.hpp"
#include <cereal/archives/binary.hpp>
#include <hibf/config.hpp>
#include <hibf/hierarchical_interleaved_bloom_filter.hpp>

void search(configuration const & config)
{
    seqan::hibf::hierarchical_interleaved_bloom_filter hibf;

    std::ifstream is{config.index_file, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(hibf);

    seqan3::sequence_file_input<seqan3::sequence_file_input_default_traits_dna> reads_file{config.reads};

    auto agent = hibf.membership_agent();
    size_t threshold = config.threshold;

    //for storing the results
    std::vector<std::string> result_for_cout;
    std::vector<std::string> result_for_txt;

    for (auto & record : reads_file)
    {
        if (record.sequence().size() < config.kmer_size)
        {
            throw std::runtime_error{"read in file is shorter than the k-mer size."};
        }

        uint8_t current_hash = determine_current_hash(config);

        auto minimiser_view =
            record.sequence()
            | seqan3::views::minimiser_hash(seqan3::ungapped{config.kmer_size}, seqan3::window_size{current_hash});

        auto & result = agent.membership_for(minimiser_view, threshold);
        agent.sort_results();

        std::string current_read = record.id() + ": [";

        for (size_t i = 0; i < result.size(); ++i)
        {
            current_read += std::to_string(result[i]);
            if (i < result.size() - 1)
                current_read += ",";
        }
        current_read += "]\n";

        // store the result in the vector
        result_for_cout.push_back(current_read);
        result_for_txt.push_back(current_read);
    }

    std::cout << "The following hits were found:\n";
    //print to console
    for (auto & record : result_for_cout)
    {
        std::cout << record;
    }
    //save to file
    std::ofstream result_out{config.search_output};
    for (auto & record : result_for_txt)
    {
        result_out << record;
    }
}
