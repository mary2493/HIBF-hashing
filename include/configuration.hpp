// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#pragma once

#include <filesystem>

struct configuration
{
    std::filesystem::path file_list_path{};
    std::filesystem::path index_output{"index"};
    uint8_t kmer_size{20u};
    bool kmer_modus{false};
    std::filesystem::path reads{};
    std::filesystem::path search_output{"output.txt"};
    std::filesystem::path index_file{};
};
