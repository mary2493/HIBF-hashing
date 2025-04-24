// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#pragma once

#include <filesystem>

enum class hash_type : uint8_t
{
    kmer,
    minimiser,
    syncmer
};

struct configuration
{
    std::filesystem::path file_list_path{};
    std::filesystem::path index_output{"index"};
    uint8_t kmer_size{20u};
    std::filesystem::path reads{};
    std::filesystem::path search_output{"output.txt"};
    std::filesystem::path index_file{};
    uint8_t error{0u};
    hash_type hash{hash_type::kmer};
    uint8_t window_size{20u};
    uint8_t s{11u};
    uint8_t t{2u};
};
