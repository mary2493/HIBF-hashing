#include "build/build.hpp"
#include <seqan3/io/sequence_file/all.hpp> 
#include <iostream>
#include <filesystem>

void build(configuration const & config)
{
    std::cout << config.file_list_path << '\n';
    std::cout << config.index_output << '\n';
    std::cout << config.kmer_size << '\n';
}
