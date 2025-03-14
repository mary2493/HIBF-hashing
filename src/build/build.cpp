#include "build/build.hpp"

#include <iostream>

void build(configuration const & config)
{
    std::cout << config.file_list_path << '\n';
    std::cout << config.index_output << '\n';
    std::cout << config.kmer_size << '\n';
}
