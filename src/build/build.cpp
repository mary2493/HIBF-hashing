#include "build/build.hpp"

#include <algorithm> // for std::all_of
#include <cctype>    // for isspace
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

#include <cereal/archives/binary.hpp>
#include <hibf/config.hpp>
#include <hibf/hierarchical_interleaved_bloom_filter.hpp>

void build(configuration const & config)
{
    std::ifstream file_list(config.file_list_path);
    std::string current_line;

    //all_bins_together: used to store each FASTA file as a separate bin
    std::vector<std::vector<uint64_t>> all_bins_together;

    //Counting the number of files in the file list
    unsigned int count_files = 0;

    //Each FASTA file is opened, and the k-mers are extracted from it.
    //These kmers are stored in all_bins_together, with each file corresponding to a "User Bin" in the HIBF
    while (std::getline(file_list, current_line))
    {

        /* Checking whether the line is empty or consists only of whitespace
        .empty() checks if the string is empty, and std::all_of checks if all characters in the string are whitespace. It returns true
        if every character is a whitespace (e.g., space, tab, newline) and false if at least one character is not a whitespace */
        if (current_line.empty() || std::all_of(current_line.begin(), current_line.end(), isspace))
        {
            throw std::runtime_error{"Empty line or invalid entry in the file list."};
        }

        seqan3::sequence_file_input<seqan3::sequence_file_input_default_traits_dna> current_fasta_file{current_line};
        std::vector<uint64_t> bin_for_kmers;

        //kmers of each file are assigned to a user bin
        for (auto & record : current_fasta_file)
        {
            //Checking if the sequence is shorter than the k-mer size
            if (record.sequence().size() < config.kmer_size)
            {
                throw std::runtime_error{"Sequence in file " + current_line + " is shorter than the k-mer size."};
            }

            //Extracting kmers from the sequence and storing them in bin_for_kmers
            auto kmers =
                record.sequence() | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{config.kmer_size}});
            bin_for_kmers.insert(bin_for_kmers.end(), kmers.begin(), kmers.end());
        }
        all_bins_together.push_back(std::move(bin_for_kmers));
        count_files++;
    }
    //The code below for building the HIBF index was copied from the https://github.com/seqan/hibf/tree/main website and adapted to the task
    //Iteration over all k-mers in the corresponding all_bins_together[user_bin_id] and insertion of these k-mers into it
    auto get_user_bin_data = [&](size_t const user_bin_id, seqan::hibf::insert_iterator it)
    {
        for (auto value : all_bins_together[user_bin_id])
            it = value;
    };

    // construct a config
    seqan::hibf::config hibf_config{.input_fn = get_user_bin_data,                   // required
                                    .number_of_user_bins = all_bins_together.size(), // required
                                    .number_of_hash_functions = 2u,
                                    .maximum_fpr = 0.05,
                                    .threads = 1u};

    // The HIBF constructor will determine a hierarchical layout for the user bins and build the filter
    seqan::hibf::hierarchical_interleaved_bloom_filter hibf{hibf_config};

    //The indices can also be stored and loaded from disk by using cereal
    std::ofstream os{config.index_output, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(hibf);

    if (count_files == 0)
    {
        throw std::runtime_error{"No valid files found in the file list."};
    }
    else
    {
        std::cout << "HIBF index built and saved to " << config.index_output << '\n';
        std::cout << "Successfully processed " << count_files << " files.";
    }
}
