// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================
//          Author: Marie Hoffmann <ozymandiaz147 AT gmail.com>
//          Manual: https://github.com/mariehoffmann/PriSeT

// Primer Settings.

// TODO: use seqan3 dna types, keep conversion table for ambiguous 1-letter encodings
#pragma once

#include <cassert>
#include <unordered_map>

//#include <seqan/basic.h>

#include "types/PrimerConfig.hpp"
#include "utilities.hpp"

namespace priset
{

// std::pair<uint64_t, uint64_t> split(TKmerID kmerID);

extern inline uint64_t complement(uint64_t const code_)
{
    auto [prefix, code] = split_kmerID(code_);
    // std::pair<uint64_t, uint64_t> prefix_code = split_kmerID(code_);
    // uint64_t prefix = prefix_code.first;
    // uint64_t code = prefix_code.second;
    uint64_t code_c = 0;
    uint8_t offset = 0;
    while (code > 1)
    {
        switch (code & 0b11)
        {
            case 0b00: code_c |= (0b11UL << (offset << 1)); break;
            case 0b01: code_c |= (0b10UL << (offset << 1)); break;
            case 0b10: code_c |= (0b01UL << (offset << 1)); break;
        }
        code >>= 2;
        ++offset;
    }
    code_c |= 1ULL << (offset << 1);
    return (prefix) ? prefix | code_c : code_c;
}

// Preserves length bits and closure.
extern inline uint64_t reverse(uint64_t const code_)
{
    auto [prefix, code] = split_kmerID(code_);
    uint64_t code_r = 0b01; // closure
    while (code > 1)
    {
        code_r <<= 2;
        code_r |= (code & 0b11);
        code >>= 2;
    }
    return prefix | code_r;
}


extern inline uint64_t reverse_complement(uint64_t const code_)
{
    auto [prefix, code] = split_kmerID(code_);
    uint64_t code_rc = 0b01; // closure
    while (code > 1)
    {
        code_rc <<= 2;
        switch (code & 0b11)
        {
            case 0b00: code_rc |= 0b11; break;
            case 0b01: code_rc |= 0b10; break;
            case 0b10: code_rc |= 0b01;
        }
        code >>= 2;
    }
    return prefix | code_rc;
}

/*
extern inline std::string reverse_complement(std::string const & seq)
{
    char rc[PRIMER_MAX_LEN];
    uint8_t i = 0;
    for (char const c : seq)
    {
        switch (c)
        {
            case 'A': rc[i] = 'T'; break;
            case 'C': rc[i] = 'G'; break;
            case 'G': rc[i] = 'C'; break;
            default: rc[i] = 'A';
        }
        ++i;
    }
    std::string seq_rc(rc, i);
    std::reverse(seq_rc.begin(), seq_rc.end());
    return seq_rc;
}
*/

    //using dna = seqan::Dna;
//!\brief DNA codes as enums.
// TODO: use seqan Dna
//enum dna {A, C, G, T, B, CGT, D, AGT, H, ACT, K, GT, M, AC, N, ACGT, R, AG, S, CG, V, ACG, W, AT, Y, CT};

/*
//!\brief Complement map of DNA codes.
std::array<dna,4> cdna = {dna::T, dna::G, dna::C, dna::A};

//!\brief String representation of DNA enums.
std::array<std::string, 26> dna2str = {"A", "C", "G", "T", "B", "CGT", "D", \
    "AGT", "H", "ACT", "K", "GT", "M", "AC", "N", "ACGT", "R", "AG", "S", "CG", \
    "V", "ACG", "W", "AT", "Y", "CT"};

//!\brief Conversion table from string representation to dna enum.
std::unordered_map<std::string, dna> str2dna = {{"A", dna::A}, {"C", dna::C}, \
      {"G", dna::G}, {"T", dna::T}, {"B", dna::B}, {"CGT", dna::CGT}, \
      {"D", dna::D}, {"AGT", dna::AGT}, {"H", dna::H}, {"ACT", dna::ACT}, \
      {"K", dna::K}, {"GT", dna::GT}, {"M", dna::M}, {"AC", dna::AC}, \
      {"N", dna::N}, {"ACGT", dna::ACGT}, {"R", dna::R}, {"AG", dna::AG}, \
      {"S", dna::S}, {"CG", dna::CG}, {"V", dna::V}, {"ACG", dna::ACG}, \
      {"W", dna::W}, {"AT", dna::AT}, {"Y", dna::Y}, {"CT", dna::CT}};

template<typename sequence_type>
static std::string dnaseq2str(sequence_type s)
{
    assert((std::is_same<dna, typename sequence_type::value_type>::value));
    std::string sstr;
    for (dna base : s) sstr.append(dna2str[base]);
    return sstr;
}

//!\brief translate into complementary string without reversing
template<typename sequence_type>
sequence_type complement(sequence_type const sequence)
{
    assert((std::is_same<dna, typename sequence_type::value_type>::value));
    sequence_type csequence = sequence;
    std::transform(csequence.begin(), csequence.end(), csequence.begin(),
                   [](dna c) -> dna { return priset::cdna[c]; });
    return csequence;
}

//!\brief bidirectional conversion table for one letter encodings.
// Attention: call always with hashed input, e.g. dna_decode[std::hash<sequence_type>("D")]
static constexpr std::array<priset::dna, 26> dna_decode
{
    [] () constexpr
    {
        std::array<priset::dna, 26> bimap{};
        bimap[priset::dna::A] = priset::dna::A;
        bimap[priset::dna::C] = priset::dna::C;
        bimap[priset::dna::G] = priset::dna::G;
        bimap[priset::dna::T] = priset::dna::T;
        bimap[priset::dna::B] = priset::dna::CGT;
        bimap[priset::dna::CGT] = priset::dna::B;
        bimap[priset::dna::D] = priset::dna::AGT;
        bimap[priset::dna::AGT] = priset::dna::D;
        bimap[priset::dna::H] = priset::dna::ACT;
        bimap[priset::dna::ACT] = priset::dna::H;
        bimap[priset::dna::K] = priset::dna::GT;
        bimap[priset::dna::GT] = priset::dna::K;
        bimap[priset::dna::M] = priset::dna::AC;
        bimap[priset::dna::AC] = priset::dna::M;
        bimap[priset::dna::N] = priset::dna::ACGT;
        bimap[priset::dna::ACGT] = priset::dna::N;
        bimap[priset::dna::R] = priset::dna::AG;
        bimap[priset::dna::AG] = priset::dna::R;
        bimap[priset::dna::S] = priset::dna::CG;
        bimap[priset::dna::CG] = priset::dna::S;
        bimap[priset::dna::V] = priset::dna::ACG;
        bimap[priset::dna::ACG] = priset::dna::V;
        bimap[priset::dna::W] = priset::dna::AT;
        bimap[priset::dna::AT] = priset::dna::W;
        bimap[priset::dna::Y] = priset::dna::CT;
        bimap[priset::dna::CT] = priset::dna::Y;
        return bimap;
    }()


};*/
}
