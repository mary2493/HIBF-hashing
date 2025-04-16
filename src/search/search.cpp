#include "search/search.hpp"

#include <iostream>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/search.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

#include "index_data.hpp"
#include <cereal/archives/binary.hpp>
#include <hibf/config.hpp>
#include <hibf/hierarchical_interleaved_bloom_filter.hpp>

void search(configuration const & config)
{
    using sequence_file_t = seqan3::sequence_file_input<seqan3::sequence_file_input_default_traits_dna>;

    myindex index{};
    index.load(config.index_file);

    sequence_file_t reads_file{config.reads};

    auto agent = index.hibf.membership_agent();
    size_t threshold = 1u;

    auto hash_adaptor =
        seqan3::views::minimiser_hash(seqan3::ungapped{index.kmer_size}, seqan3::window_size{index.window_size});

    //for storing the results
    std::vector<std::string> results;
    std::string current_read{};

    for (auto & record : reads_file)
    {
        if (record.sequence().size() < index.window_size)
        {
            throw std::runtime_error{"read in file is shorter than the k-mer/window size."};
        }

        auto minimiser_view = record.sequence() | hash_adaptor;

        auto & result = agent.membership_for(minimiser_view, threshold);
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
