// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#pragma once

#include "configuration.hpp"

inline uint8_t determine_current_hash(configuration const & config)
{
    if (config.hash == hash_type::kmer)
        return config.kmer_size;
    else if (config.hash == hash_type::minimiser)
        return config.window_size;
    else
        throw std::runtime_error{"Syncmer support is not yet implemented. Please use kmer or minimiser."};
}