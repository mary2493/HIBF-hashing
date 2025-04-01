#include "search/search.hpp"
#include <iostream>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/search.hpp>

using namespace seqan3::literals;

void search(configuration const &)
{
   
        std::ifstream is{config.index_file, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(hibf);
     
        seqan3::debug_stream << "=====   Running on a single text   =====\n"
                             << "The following hits were found:\n";
     
        for (auto && result : search("GCT"_dna4, index))
            seqan3::debug_stream << result << '\n';
    
     
    void run_text_collection()
    {
        std::vector<seqan3::dna4_vector> text{"CGCTGTCTGAAGGATGAGTGTCAGCCAGTGTA"_dna4,
                                              "ACCCGATGAGCTACCCAGTAGTCGAACTG"_dna4,
                                              "GGCCAGACAACCCGGCGCTAATGCACTCA"_dna4};
        seqan3::fm_index index{text};
     
        seqan3::debug_stream << "===== Running on a text collection =====\n"
                             << "The following hits were found:\n";
     
        for (auto && result : search("GCT"_dna4, index))
            seqan3::debug_stream << result << '\n';
    }
     



}
