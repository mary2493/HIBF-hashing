// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#pragma once

#include "configuration.hpp"
#include <cereal/archives/binary.hpp>
#include <hibf/config.hpp>
#include <hibf/hierarchical_interleaved_bloom_filter.hpp>

class myindex
{
public:
    uint8_t kmer_size{};
    uint8_t window_size{};
    uint8_t s{};
    uint8_t t{};
    hash_type hash{};
    seqan::hibf::hierarchical_interleaved_bloom_filter hibf{};

    myindex() = default;
    myindex & operator=(myindex const &) = default;
    myindex(myindex const &) = default;
    myindex(myindex &&) = default;
    myindex & operator=(myindex &&) = default;
    ~myindex() = default;

    explicit myindex(configuration const & config, seqan::hibf::hierarchical_interleaved_bloom_filter index) :
        kmer_size{config.kmer_size},
        window_size{config.window_size},
        s{config.s},
        t{config.t},
        hash{config.hash},
        hibf{std::move(index)}
    {}

    void store(std::filesystem::path const & path) const
    {
        std::ofstream fout{path};
        cereal::BinaryOutputArchive oarchive{fout};
        oarchive(*this);
    }

    void load(std::filesystem::path const & path)
    {
        std::ifstream fin{path};
        cereal::BinaryInputArchive iarchive{fin};
        iarchive(*this);
    }

    template <typename archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        archive(kmer_size);
        archive(window_size);
        archive(s);
        archive(t);
        archive(hash);
        archive(hibf);
    }
};
