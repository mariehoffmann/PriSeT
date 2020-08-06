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
#include "common.hpp"
#include "types/GenMapTypes.hpp"
#include "types/PrimerConfig.hpp"
#include "utilities.hpp"

namespace priset
{

uint64_t get_code(uint64_t const code_, uint64_t mask);

/* Encode a single sequence as a 64 bit integer.
 * Details: encoding schme is \sum_i 4^i*x_i, starting with the first character
 * (little endian) and x_i being the 2 bit representation of 'A' (=0), 'C' (=1),
 * 'G' (=2), and 'T' (=4) ..., 3 = 'G'. E.g., ACGT is encoded as 0*4^0 + 1*4^1 + 2*4^2 + 3*4^2.
 * Non-zero encoded character ('C') is added because of flexible sequence lengths
 * and therefore the necessity to differentiate 'XA' from 'XAA'.
 * Since multiple kmers may start at one position in the reference (which means that they all
 * share the same prefix), the code stores in its 12 highest bits flags for which length
 * the code represents. E.g. if the 1st bit is set, the code represents a string
 * with the length of the shortest possible primer length, if the 5th bit is set, the code also represents
 * a string of minimal primer length plus 4, and so forth.
 *
 * Bit layout:
 * bits [0 .. |seq|-1]  complete sequence
 * bit [|seq|]          closure symbol 'C'
 * bits [60:64]         lower sequence length bound in case of variable length
 */

uint64_t dna_encoder(seqan::String<dna> const & seq)
{
    uint64_t code(0);
     for (uint64_t i = 0; i < seqan::length(seq); ++i)
     {
         switch (char(seqan::getValue(seq, i)))
         {
             case 'C': code |= 1ULL; break;
             case 'G': code |= 2ULL; break;
             case 'T': code |= 3ULL;
         }
         code <<= 2;

     }
     code >>= 2;
     code |= 1ULL << uint64_t(seqan::length(seq) << 1); // stop symbol 'C' = 1
     return code;
}

// with length bit in prefix
uint64_t dna_encoder_with_lbit(seqan::String<priset::dna> const & seq)
{
    uint64_t code(0);
    auto l = seqan::length(seq);
    for (uint64_t i = 0; i < l; ++i)
    {
         switch (char(seqan::getValue(seq, i)))
         {
             case 'C': code |= 1ULL; break;
             case 'G': code |= 2ULL; break;
             case 'T': code |= 3ULL;
         }
         code <<= 2;
     }
     // revert last shift
     code >>= 2;
     // add length bit and stop symbol 'C' = 1
     code |= 1ULL << (63 - (l - KAPPA_MIN)) | (1ULL << uint64_t(l << 1));
     return code;
}

// Return dna sequence given code and length mask.
std::string dna_decoder(uint64_t const _code, uint64_t const mask)
{
    if (_code == 0ULL)
        throw std::invalid_argument("ERROR: invalid argument for decoder(0).");
    uint64_t code = _code;
    code = get_code(_code, mask);
    std::array<char, 4> alphabet = {'A', 'C', 'G', 'T'};
    uint8_t const n = (63 - __builtin_clzl(code)) >> 1;
    char * seq_ptr = new char[n];
    for (uint8_t i = 1; i <= n; ++i, code >>= 2)
        seq_ptr[n - i] = alphabet[3 & code];
    delete[] seq_ptr;
    return std::string(seq_ptr, n);
}

// return full length sequence, ignore variable length info in leading bits if present.
// Note: kmer code is only length trimmed when a mask is given.
std::string dna_decoder(uint64_t const _code)
{
    if (_code == 0ULL)
        throw std::invalid_argument("ERROR: invalid argument for decoder(0).");
    uint64_t code = _code;
    code &= ~PREFIX_SELECTOR;
    std::array<char, 4> alphabet = {'A', 'C', 'G', 'T'};
    uint8_t const n = (63 - __builtin_clzl(code)) >> 1;
    char * seq_ptr = new char[n];
    for (uint8_t i = 1; i <= n; ++i, code >>= 2)
        seq_ptr[n - i] = alphabet[3 & code];
    delete[] seq_ptr;
    return std::string(seq_ptr, n);
}

extern inline uint64_t complement(uint64_t const _code)
{
    auto [prefix, code] = split_kmerID(_code);
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
extern inline uint64_t reverse(uint64_t const _code)
{
    auto [prefix, code] = split_kmerID(_code);
    uint64_t code_r = 0b01; // closure
    while (code > 1)
    {
        code_r <<= 2;
        code_r |= (code & 0b11);
        code >>= 2;
    }
    return prefix | code_r;
}


extern inline uint64_t reverse_complement(uint64_t const _code)
{
    auto [prefix, code] = split_kmerID(_code);
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
