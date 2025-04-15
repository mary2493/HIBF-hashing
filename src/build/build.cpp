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

#include <cereal/archives/binary.hpp>
#include <hibf/config.hpp>
#include <hibf/hierarchical_interleaved_bloom_filter.hpp>

void build(configuration const & config)
{
    using sequence_file_t = seqan3::sequence_file_input<seqan3::sequence_file_input_default_traits_dna>;

    std::ifstream file_list(config.file_list_path);
    std::string current_line;

    std::vector<std::string> user_bin_paths;

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

    auto get_user_bin_data = [&](size_t const user_bin_id, seqan::hibf::insert_iterator it)
    {   
        sequence_file_t fin{user_bin_paths[user_bin_id]};
        for (auto & record : fin)
        {
            if (record.sequence().size() < config.kmer_size)
                throw std::runtime_error{"Sequence in " + user_bin_paths[user_bin_id]
                                         + " is shorter than the k-mer size."};

            std::ranges::copy(record.sequence()
                                  | seqan3::views::minimiser_hash(seqan3::ungapped{config.kmer_size},
                                                                  seqan3::window_size{config.kmer_size}),
                              it);
        }
    };
    // construct a config
    seqan::hibf::config hibf_config{.input_fn = get_user_bin_data,                // required
                                    .number_of_user_bins = user_bin_paths.size(), // required
                                    .number_of_hash_functions = 2u,
                                    .maximum_fpr = 0.05,
                                    .threads = 1u};

    // The HIBF constructor will determine a hierarchical layout for the user bins and build the filter
    seqan::hibf::hierarchical_interleaved_bloom_filter hibf{hibf_config};

    //The indices can also be stored and loaded from disk by using cereal
    std::ofstream os{config.index_output, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(hibf);

    std::cout << "HIBF index built and saved to " << config.index_output << "\n";
    std::cout << "Successfully processed " << user_bin_paths.size() << " files.\n";
}
