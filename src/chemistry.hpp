// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================
//          Author: Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
//          Manual: https://github.com/mariehoffmann/PriSeT

// Global Functions for Primer Constraint Checking.

#pragma once

#include <algorithm>
#include <bitset>
#include <cassert>
#include <cmath>
#include <iostream>
#include <functional>
#include <vector>

#include <seqan/basic.h>

#include "dna.hpp"
#include "primer_cfg_type.hpp"
#include "types.hpp"
#include "utilities.hpp"

namespace priset  //::chemistry TODO: introduce chemistry namespace
{

std::string dna_decoder(uint64_t code, uint64_t const mask);
// std::pair<uint64_t, uint64_t> split(TKmerID kmerID);
extern inline uint64_t complement(uint64_t const code);
extern inline uint64_t reverse(uint64_t const code);
extern inline uint64_t reverse_complement(uint64_t const code);
extern inline void trim_to_true_length(TKmerID &);

// Difference in melting temperatures (degree Celsius) according to Wallace rule.
extern inline float dTm(TKmerID const kmerID1, TKmerID const mask1, TKmerID const kmerID2, TKmerID const mask2)
{
    if (!(kmerID1 & mask1) || !(kmerID2 & mask2))
        std::cerr << "Error: target mask bit not set\n";

    // remove length mask
    TKmerID code1 = kmerID1 & ~PREFIX_SELECTOR;
    TKmerID code2 = kmerID2 & ~PREFIX_SELECTOR;

    uint8_t enc_l1 = WORD_SIZE - 1 - __builtin_clzl(code1);  // encoded length 2bit
    uint8_t mask_l1 = (__builtin_clzl(mask1) + PRIMER_MIN_LEN) << 1;   // selected length 2bit

    uint8_t enc_l2 = WORD_SIZE - 1 - __builtin_clzl(code2); // encoded length
    uint8_t mask_l2 = (__builtin_clzl(mask2) + PRIMER_MIN_LEN) << 1;      // selected length

    if (mask_l1 > enc_l1 || mask_l2 > enc_l2)
        std::cerr << "ERROR: largest encoded kmer length undershoots target length!\n";

    // trim kmer if exceeding encoded length
    code1 >>= enc_l1 - mask_l1;
    code2 >>= enc_l2 - mask_l2;

    int8_t AT{0};
    int8_t CG{0};

    if (!(__builtin_clzl(code1) % 2))
        std::cerr << "ERROR: expected odd number of leading zeros\n";

    if (!(__builtin_clzl(code2) % 2))
        std::cerr << "ERROR: expected odd number of leading zeros\n";

    while (code1 != 1)
    {
        (!(code1 & 3) || (code1 & 3) == 3) ? ++AT : ++CG;
        code1 >>= 2;
    }
    while (code2 != 1)
    {
        (!(code2 & 3) || (code2 & 3) == 3) ? --AT : --CG;
        code2 >>= 2;
    }
    return std::abs(2 * AT + 4 * CG);
}

// AT   AT = 2
//
// CG      CG = -2

/*
 * Compute melting temperature of an encoded oligomer. The oligomer is trimmed
 * to the length indicated by bit mask, otherwise (mask = 0) to the smallest encoded
 * length.
 * \Details
 * Compute the melting temperature applying the simple Wallace rule:
 *                            Tm = 2AT + 4CG
 * The difference in temperature estimates between the Wallace and the nearest
 * neighbour method (considered as more precise) rarely exceeds 4.5 K. This
 * maximum estimatation error needs to be considered when filtering and analyzing
 * potential primer pairs.
 */
extern inline uint8_t Tm(TKmerID const kmerID, uint64_t const mask)
{
    auto [prefix, code] = split_kmerID(kmerID);
    auto target_l = PRIMER_MIN_LEN;
    target_l += (prefix & !mask) ? __builtin_clzll(prefix) : __builtin_clzll(mask);
    auto enc_l = (WORD_SIZE - 1 - __builtin_clzll(code)) >> 1;
    code >>= (enc_l - target_l) << 1;
    uint8_t AT = 0;
    uint8_t CG = 0;
    while (code != 1)
    {
        switch (code & 3)
        {
            case 0:
            case 3: ++AT; break;
            case 1:
            case 2: ++CG;
        }
        code >>= 2;
    }
    return (AT << 1) + (CG << 2);
}

// CG content computation for output
extern inline float CG(TKmerID const kmerID, uint64_t const mask)
{
    auto [prefix, code] = split_kmerID(kmerID);
    auto target_l = PRIMER_MIN_LEN;
    target_l += (prefix & !mask) ? __builtin_clzll(prefix) : __builtin_clzll(mask);
    auto enc_l = (WORD_SIZE - 1 - __builtin_clzll(code)) >> 1;
    code >>= (enc_l - target_l) << 1;
    enc_l = (WORD_SIZE - 1 - __builtin_clzll(code)) >> 1;
    uint8_t CG = 0;
    while (code != 1)
    {
        switch (code & 3)
        {
            case 1:
            case 2: ++CG;
        }
        code >>= 2;
    }
    return float(CG)/float(enc_l);;
}

//!\brief Salt-adjusted method to compute melting temperature of primer sequence.
// input primer:string sequence, Na:float molar Natrium ion concentration
extern inline float primer_melt_salt(TKmerID code, float const Na)
{
    uint8_t CG = 0;
    uint8_t seq_len = 0;
    std::array<char, 4> decodes = {'A', 'C', 'G', 'T'};
    while (code != 1)
    {
        switch(decodes[3 & code])
        {
            case 'C':
            case 'G': ++CG;
        }
        code >>= 2;
        ++seq_len;
    }
    return 100.5 + 41.0 * CG / seq_len - 820.0 / seq_len + 16.6 * std::log10(Na);
}

/* !\brief Check for low energy secondary structures.
 * Returns true if stability of a hairpin structure above the annealing temperature is improbable.
 * It is assumed that the annealing temperature
 * is 5°C below the Tｍ of the sequence.
 * TODO: implement Zuker Recursion
 */
/*bool filter_hairpin(primer_cfg_type const & primer_cfg, seqan::String<priset::dna> const & sequence)
{
    float Ta = get_Tm(primer_cfg, sequence) - 5;
    // check for hairpins

        // check its melting temperature
    return true;
}*/

uint64_t get_code(uint64_t const code_, uint64_t mask);

template<typename uint_type>
std::string bits2str(uint_type i);

std::string kmerID2str(TKmerID kmerID);

/*
 * This filter is one of many filters applied sequentially. Therefore, it is
 * possible that the prefix is reset completely by previous operations.
 * TODO: shorten di-nucleotide runs to at most 4 due to self-annealing
 * Filters also a subset of self-annealing structures.
 */
extern inline void filter_repeats_runs(TKmerID & kmerID)
{
    // std::cout << "\tenter filter_repeats_runs ...\n";
    auto [prefix, code] = split_kmerID(kmerID);
    if (!code)
        throw std::invalid_argument("Expected kmerID non zero!");
    if (!prefix)
    {
        // std::cout << "prefix is already 0 \n";
        return;
    }
    uint64_t const tail_selector_10 = (1 << 10) - 1;
    uint64_t const tail_selector_20 = (1 << 20) - 1;
    kmerID = code; // save trimmed kmer-code part
    uint64_t const k = (63 - __builtin_clzl(code)) >> 1;
    for (uint64_t i = 0; i < k - 4; ++i) // up-to four, i.e. 8 bits consecutive characters permitted
    {
        uint64_t tail_10 = tail_selector_10 & code;
        // XOR tail with A_4, C_4, G_4, T_4
        if ((tail_10 == 0b00000000) || (tail_10 == 0b01010101) || (tail_10 == 0b10101010) || (tail_10 == 0b11111111))
        {
            // note that (x >> 64) won't be executed
            uint64_t offset = 64 - (std::max(PRIMER_MIN_LEN, k - i) - PRIMER_MIN_LEN);
            // delete all k bits ≥ len_selector
            prefix = (offset == 64) ? 0 : (prefix >> offset) << offset;
            // std::cout << "\t4-run found\n";
        }
        // at least 10 characters left to test for di-nucleotide repeats of length 5
        if (k - i > 9)
        {
            uint64_t tail_20 = code & tail_selector_20;
            // AT_5, TA_5, AC_5, CA_5, AG_5, GA_5, CG_5, GC_5, CT_5, TC_5, GT_5, TG_5
            if ( (tail_20 == 0b00110011001100110011) || (tail_20 == 0b11001100110011001100) ||
                 (tail_20 == 0b00010001000100010001) || (tail_20 == 0b01000100010001000100) ||
                 (tail_20 == 0b00100010001000100010) || (tail_20 == 0b10001000100010001000) ||
                 (tail_20 == 0b01100110011001100110) || (tail_20 == 0b10011001100110011001) ||
                 (tail_20 == 0b01110111011101110111) || (tail_20 == 0b11011101110111011101) ||
                 (tail_20 == 0b10111011101110111011) || (tail_20 == 0b11101110111011101110))
            {
                // delete all k bits ≥ len_selector
                uint64_t offset = 64 - (std::max(PRIMER_MIN_LEN, k - i) - PRIMER_MIN_LEN);
                // delete all k bits ≥ len_selector
                prefix = (offset == 64) ? 0 : (prefix >> offset) << offset;
                // std::cout << "\tdi-nucleotide run found\n";
            }
        }
        if (!prefix)
        {
            // std::cout << "\tprefix is 0\n";
            break;
        }

        code >>= 2;
    }
    // add updated prefix

    kmerID |= prefix;
    // std::cout << "\tored kmerID = " << kmerID2str(kmerID) << std::endl;
}

// Check if not more than 3 out of the 5 last bases at the 3' end are CG.
// DNA sense/'+': 5' to 3' tests last 5 bps, antisense/'-': 3' to 5' tests first 5 bps
// Returns false if constraint is violated.
extern inline bool filter_CG_clamp(TKmerID const kmerID, char const sense, uint64_t const mask = 0)
{
    auto [prefix, code] = split_kmerID(kmerID);
    uint64_t enc_l = WORD_SIZE - __builtin_clzll(code) - 1; // in bits
    if (sense == '-')
        code >>= enc_l - 10; // shift right to have prefix of length 10
    else
    {
        if (mask)
        {
            uint8_t target_l = (PRIMER_MIN_LEN + __builtin_clzl(mask)) << 1; // in bits
            code >>= enc_l - target_l; // delete prefix and trim to target length
        }
    }
    return  ((((code & 3) | ((code & 3) + 1)) == 3) +
            ((((code >> 2) & 3) | (((code >> 2) & 3) + 1)) == 3) +
            ((((code >> 4) & 3) | (((code >> 4) & 3) + 1)) == 3) +
            ((((code >> 6) & 3) | (((code >> 6) & 3) + 1)) == 3) +
            ((((code >> 8) & 3) | (((code >> 8) & 3) + 1)) == 3) <= 3) ? true : false;
}

// https://www.linuxquestions.org/questions/programming-9/can-i-predefine-constant-array-in-c-823350/
/*
 * Bit map for (A|T)^3 patterns:
 * 0b000000 = 0b000011 = 0b001100 = 0b001111 = 0b110000 = 0b110011 = 0b111100 = 0b111111 = 1
 */
inline uint8_t const* WWW_tails()
{

    static uint8_t const tails[64] = {  1, 0, 0, 1, 0, 0, 0, 0,    // 0
                                        0, 0, 0, 0, 1, 0, 0, 1,    // 8
                                        0, 0, 0, 0, 0, 0, 0, 0,    // 16
                                        0, 0, 0, 0, 0, 0, 0, 0,    // 24
                                        0, 0, 0, 0, 0, 0, 0, 0,    // 32
                                        0, 0, 0, 0, 0, 0, 0, 0,    // 40
                                        1, 0, 0, 1, 0, 0, 0, 0,    // 48
                                        0, 0, 0, 0, 1, 0, 0, 1};    // 56

    /* not implemented yet in gnu gcc compiler mac distribution: designation
    static uint8_t const tails[64] =
    {
        [0b000000] = 1, [0b000011] = 1, [0b001100] = 1, [0b001111] = 1,
        [0b110000] = 1, [0b110011] = 1, [0b111100] = 1, [0b111111] = 1
    };
    */
    return tails;
}

/*
 * Filter for tails with a (A|T)^3 pattern.
 */
extern inline bool filter_WWW_tail(TKmerID const kmerID, char const sense,
    uint64_t const mask = 0)
{
    if (sense != '+' && sense != '-')
        throw std::invalid_argument("Expected sense '-' or '+'!");
    auto [prefix, code] = split_kmerID(kmerID);
    uint64_t encoded_len = WORD_SIZE - __builtin_clzll(code) - 1;
    if (sense == '-') // reverse complement!
    {
        code >>= encoded_len - 6; // shift right to have prefix of length 10
        return !WWW_tails()[code & 0b111111];
    }
    if (mask)
    {
        uint64_t target_len = (PRIMER_MIN_LEN + __builtin_clzl(mask)) << 1;
        code >>= encoded_len - target_len; // delete prefix and trim to target length
    }
    return !WWW_tails()[code & 0b111111];
}

extern inline void filter_annealing_disconnected(TKmerID &, TKmerID &);

extern inline void delete_length_bits(TKmerID & kmerID, uint8_t l);

// l1, l2 refer to kmer lengths including dimerizing 4-mer
extern inline bool annealing_helper(TKmerID & kmerID1, uint8_t l1, TKmerID & kmerID2, uint8_t l2)
{

    // std::cout << "l1 = " << int(l1) << ", l2 = " << int(l2) << ", kmerID2 ? " << (kmerID2 ? 1 : 0)  << std::endl;
    if (!kmerID2) // first 4-mer occurrence significant
        delete_length_bits(kmerID1, std::min(l1, l2));
    else
    {
        // case: l1 <= 16 and l2 > 16 or l1 > 16 and l2 > 16
        if (l2 > (PRIMER_MIN_LEN << 1) && l2 > l1)
        {
            // std::cout << "l1 > primer_min_len\n";
            delete_length_bits(kmerID2, l2);
        }
        // case: l1 > 16 and l2 <= 16 or l1 > 16 and l2 > 16
        if (l1 > (PRIMER_MIN_LEN << 1) && l1 > l2)
        {
            // std::cout << "l1 > primer_min_len\n";
            delete_length_bits(kmerID1, l1);
        }
        // case: l1 <= 16 and l2 <= 16
        if (l1 <= (PRIMER_MIN_LEN << 1) && l2 <= (PRIMER_MIN_LEN << 1))
        {
            // std::cout << "both kmers in prefixes\n";
            delete_length_bits(kmerID1, l1);
            delete_length_bits(kmerID2, l2);
        }
    }
    return (!(kmerID1 & PREFIX_SELECTOR) || !(kmerID2 & PREFIX_SELECTOR)) ? 1 : 0;
}

/*
 * Self- and cross-dimerization test (if 2nd kmer is not zero). 4-mers ending before
 * position PRIMER_MIN_LEN result in deletion of length mask, otherwise only
 * affected lengths are erased. In the second part non-consecutive annealing is
 * tested - more than 8 annealing positions are considered to be critical.
 */
 // TODO: cross-annealing 4-mer in suffix 17,17 only excludes 17-17, but still allows for
 // combining 16-17 and vice versa, and 16 - 16, deleting length bit 17 in one kmer may exclude
 // one ore more possible combinations, better directly receeive combiner struct to set
 // combination bits instead of post-processing length bits of both k-mers
extern inline void filter_annealing_connected(TKmerID & kmerID1, TKmerID & kmerID2 = NULL_TKMERID)
{

    uint64_t const fourmer_mask = 0b11111111;
    auto [prefix1, code1] = split_kmerID(kmerID1);
    auto [prefix2, code2] = (kmerID2) ? split_kmerID(kmerID2) : split_kmerID(kmerID1);
    if (!prefix1 || !prefix2)
        return;
    trim_to_true_length(kmerID1);
    trim_to_true_length(kmerID2);
    uint64_t const code2_rev = reverse(code2);
    uint8_t const l1 = WORD_SIZE - __builtin_clzll(code1) - 1;
    uint8_t const l2 = WORD_SIZE - __builtin_clzll(code2) - 1;
    uint64_t prefix;
    uint8_t overlap;
    uint64_t overlap_mask;
    int ctr = 32;
    // test overlap positions [0, 0:l1-2] @code1 against [l2-8:max(0,l2-l1+2), l2] @code2
    for (overlap = 8; (overlap_mask = (1ULL << overlap) - 1) <= (1ULL << (l1 - 2)) - 1; ++++overlap) // TODO: test invariant
    {
        if (ctr < 0)
            exit(0);
        --ctr;
        // overlap_mask = (1ULL << overlap) - 1;
        // prefix of code1 and suffix of code2
        prefix = (code1 >> (l1 - overlap));
        // suffix = code2; // & overlap_mask; // needed?
        // suffix_rev = code2_rev;// & overlap_mask;
        // handle case l1 > l2, trim leading bits from 1st code
        if (kmerID2 && overlap > l2 && l2 < l1)
            overlap_mask = (1 << l2) - 1;
        // std::cout << "overlap_mask = " << std::bitset<36>(overlap_mask).to_string() << std::endl;
        auto x = (prefix ^ code2) & overlap_mask;
        auto y = (prefix ^ code2_rev) & overlap_mask;

        uint8_t i = 0;
        while (x >= fourmer_mask || y >= fourmer_mask)
        {
            if ((x & fourmer_mask) == fourmer_mask)
            {

                if (annealing_helper(kmerID1, overlap - i, kmerID2, l2 - i))
                    return;
            }
            if ((y & fourmer_mask) == fourmer_mask)
            {
                if (annealing_helper(kmerID1, overlap - i, kmerID2, i + 8))
                    return;
            }
            x >>= 2;
            y >>= 2;
            ++++i;
        }
    }
    // test for cross- and inverse-dimerization, i.e. overlap positions
    // [0:l1-8, l1] @code against [l2-l1:0, l2:8]
    if (kmerID1 & PREFIX_SELECTOR)
    {
        overlap = (l1 < l2) ? l1 : l2;
        uint8_t code2_prefix_ctr = (l1 < l2) ? l2 - l1 : 0;
        uint8_t offset = 0;
        while ((overlap_mask = (1ULL << overlap) - 1) >= fourmer_mask)
        {
            if (kmerID2)
            {
                auto x = (code1 ^ (code2 >> offset)) & overlap_mask;
                uint8_t i = 0;
                while (x >= fourmer_mask)
                {
                    if ((x & fourmer_mask) == fourmer_mask)
                    {
                        if (annealing_helper(kmerID1, l1 - i, kmerID2, l2 - offset - i))
                            return;
                    }
                    x >>= 2;
                    ++++i;
                }
            }
            // cross- or self-annealing with reversed sequence
            auto y = (code1 ^ (code2_rev >> offset)) & overlap_mask;
            uint8_t i = 0;
            while (y >= fourmer_mask)
            {
                if ((y & fourmer_mask) == fourmer_mask)
                {
                    if (annealing_helper(kmerID1, l1 - i, kmerID2, offset + i + 8))
                        return;
                }
                y >>= 2;
                ++++i;
            }
            // process code2 prefix if code2 longer than code1
            if (l1 < l2 && code2_prefix_ctr)
                ----code2_prefix_ctr;
            else // continue shortening overlap
                ----overlap;
            ++++offset;
        }
    }
}

// check for disconnected self-annealing positions <= 8
// code1 <--------------
//         |||||||||||||
// code2   -------------->
//         |||||||||||||
// code3 >--------------
extern inline void filter_annealing_disconnected(TKmerID & kmerID1, TKmerID & kmerID2 = NULL_TKMERID)
{
    if (!(kmerID1 & PREFIX_SELECTOR) && !(kmerID2 & PREFIX_SELECTOR))
        return;
    // std::cout << "l1 before trim: " << WORD_SIZE - __builtin_clzll(kmerID1&~PREFIX_SELECTOR) - 1 << std::endl;

    trim_to_true_length(kmerID1);
    // std::cout << "l1 after trim: " << WORD_SIZE - __builtin_clzll(kmerID1&~PREFIX_SELECTOR) - 1 << std::endl;
    if (kmerID2)
        trim_to_true_length(kmerID2);
    auto [prefix1, code1] = split_kmerID(kmerID1);
    auto [prefix2, code2] = (kmerID2) ? split_kmerID(kmerID2) : split_kmerID(kmerID1);
    uint64_t code2_rev = reverse(code2);
    uint8_t l1 =  WORD_SIZE - __builtin_clzll(code1) - 1;
    uint8_t l2 =  WORD_SIZE - __builtin_clzll(code2) - 1;

    // std::cout << "l1 = " << int(l1) << ", l2 = "  << int(l2) << std::endl;
    // std::cout << std::bitset<64>(code1) << std::endl;
    uint8_t overlap = 20;
    uint64_t overlap_mask;
    uint64_t prefix, x, y;

    // Due to symmetry we need to test only the first half of possible overlapping positions
    // std::cout << std::bitset<50>((1ULL << (std::min(l1, l2) - 2)) - 1) << std::endl;
    // std::cout << "init overlap mask = " << std::bitset<50>((1ULL << (overlap << 1)) - 1) << std::endl;
    // std::cout << "maximal forward overlap = " << std::min(l1, l2) - 2 << std::endl;
    // std::cout << "initial overlap = " << overlap << std::endl;
    for (overlap = 20; (overlap_mask = (1ULL << overlap) - 1) <= (1ULL << (std::min(l1, l2) - 2)) - 1; ++++overlap)
    {
        // prefix of code1 and suffix of code2
        prefix = (code1 >> (l1 - overlap));
        // creates 11 patterns with even offsets for complementary bases
        x = (prefix ^ code2) & overlap_mask;
        // std::cout << "x1 = " << std::bitset<50>(x) << std::endl;

        y = (prefix ^ code2_rev) & overlap_mask;

        // count ones every two bits and store in 2 bits, 011011 -> 010110
        x = x - ((x >> 1) & 0x1555555555555);
        // std::cout << "x2 = " << std::bitset<50>(x) << std::endl;

        // filter every second set bit
        x &= 0xAAAAAAAAAAAA & overlap_mask; // 0b10_{50}
        // std::cout << "x3 = " << std::bitset<50>(x) << std::endl;
        auto annealings_x = __builtin_popcountll(x);
        // std::cout << "overlap = " << int(overlap) << ", annealings_x = " << annealings_x << std::endl;
        if (annealings_x >= 8)
            annealing_helper(kmerID1, overlap, kmerID2, l2);
        y = y - ((y >> 1) & 0x1555555555555);
        // filter every second set bit
        y &= 0xAAAAAAAAAAAA; // 0b10_{42}
        auto annealings_y = __builtin_popcountll(y);
        if (annealings_y >= 8)
            annealing_helper(kmerID1, overlap, kmerID2, l2);
    }

    // std::cout << "-------------------------------------------------------\n";
    if (kmerID1 & PREFIX_SELECTOR)
    {
        // std::cout << "l1 = " << WORD_SIZE - __builtin_clzll(code1) - 1 << std::endl;
        overlap = (l1 <= l2) ? l1 : l2;
        // std::cout << "initial overlap = " << int(overlap) << std::endl;
        while (overlap >= 20)
        {
            overlap_mask = (1ULL << overlap) - 1;
            if (kmerID2) // cross-annealing with code2 forward
            {
                x = (code1 ^ code2) & overlap_mask;
                x = x - ((x >> 1) & 0x1555555555555);
                x &= 0xAAAAAAAAAAAA & overlap_mask; // 0b10_{50}
                // std::cout << "x3 = " << std::bitset<50>(x) << std::endl;
                auto annealings_x = __builtin_popcountll(x);
                if (annealings_x >= 8)
                {
                    kmerID1 &= ~PREFIX_SELECTOR;
                    kmerID2 &= ~PREFIX_SELECTOR;
                    break;
                }
            }
            // std::cout << "y = " << std::bitset<50>(y) << std::endl;
            y = (code1 ^ code2_rev) & overlap_mask;
            // std::cout << "y = " << std::bitset<50>(y) << std::endl;
            y = y - ((y >> 1) & 0x1555555555555);
            // std::cout << "y = " << std::bitset<50>(y) << std::endl;
            y &= 0xAAAAAAAAAAAA & overlap_mask; // 0b10_{42}
            // std::cout << "y = " << std::bitset<50>(y) << std::endl;
            auto annealings_y = __builtin_popcountll(y);
            // std::cout << "overlap = " << int(overlap) << ", annealings = " << int(annealings_y) << std::endl;
            // std::cout << "reverse 2nd half: " << std::bitset<50>(y) << std::endl;
            // TODO: in case l1 < l2, for (l2-l1) iterations l2 is smaller, but then increases to l2
            if (annealings_y >= 8)
            {
                kmerID1 &= ~PREFIX_SELECTOR;
                kmerID2 &= ~PREFIX_SELECTOR;
                break;
            }
            code2 >>= 2;
            // std::cout << "code2 = " <<  kmerID2str(code2) << std::endl;
            overlap = std::min(WORD_SIZE - __builtin_clzll(code2) - 1, int(l1));
        }
    }
}

/*
extern inline void filter_annealing(TKmerID & kmerID) // function overload for self-annealing
{
    static TKmerID const kmerID2{0};
    filter_annealing(kmerID, kmerID2);
}
*/

// test forward and reverse cross annealing between two kmers. Two kmers cross-anneal
// if there are at least 4 nts complementary to each other.
extern inline void filter_cross_annealing(TKmerID & kmerID1, TKmerID & kmerID2)
{
    auto [prefix1, code1] = split_kmerID(kmerID1);
    auto [prefix2, code2] = split_kmerID(kmerID2);
    if (!prefix1 || !prefix2)
        return;
    // build 4-mer position map of kmerID1
    uint8_t infix_mask = 0b11111111;
    uint8_t four[512];
    std::memset(four, std::numeric_limits<uint8_t>::max(), 512);
    while (code1 > (1 << 8))
    {
        four[infix_mask & code1] = uint8_t(WORD_SIZE - __builtin_clzl(code1) - 1);
        code1 >>= 2;
    }
    // forward/reverse annealing of kmerID2
    while (code2 > (1 << 8))
    {
        // same sense homology
        uint64_t infix = (~code2 & infix_mask);
        // antisense homology
        uint64_t infix_rev = ((infix & 0b11) << 6) + ((infix & 0b1100) << 2) + ((infix & 0b110000) >> 2) + ((infix & 0b11000000) >> 6);
        if (four[infix] != std::numeric_limits<uint8_t>::max() || four[infix_rev] != std::numeric_limits<uint8_t>::max())
        {
            // select earliest annealing pattern
            uint8_t l1 = std::min(four[infix], four[infix_rev]) + 8;
            uint8_t l2 = WORD_SIZE - __builtin_clzl(code2) - 1;
            // annealing in both prefix makes them incombinable
            if (l1 <= PRIMER_MIN_LEN && l2 <= PRIMER_MIN_LEN)
            {
                kmerID1 &= ~PREFIX_SELECTOR;
                kmerID2 &= ~PREFIX_SELECTOR;
                return;
            }
            else
            {
                // kmerID1 remains combinable with prefix of kmerID2 or l2 prefix needs to be trimmed
                if (l1 <= PRIMER_MIN_LEN || l2 > PRIMER_MIN_LEN)
                    kmerID2 = (kmerID2 & ~PREFIX_SELECTOR) + (prefix2 & ~((ONE_LSHIFT_63 >> (l2 - PRIMER_MIN_LEN - 1)) - 1));
                if (l2 <= PRIMER_MIN_LEN || l1 > PRIMER_MIN_LEN) // kmerID2 remains combinable with prefix of kmerID1
                    kmerID1 = (kmerID1 & ~PREFIX_SELECTOR) + (prefix1 & ~((ONE_LSHIFT_63 >> (l1 - PRIMER_MIN_LEN - 1)) - 1));
            }
        }
        code2 >>= 2;
    }
}

/*
 * Filter kmers based on their chemical properties regardless of their pairing.
 * This is a metafunction performing several single pass checks:
 *   1. melting tempaerature
 *   2. CG content
 *   3. di-nucleotide repeats
 *   4. consecutive runs
 * Length bits in prefix of kmerID are reset in case of not passing.
 * This function only modifies the prefix.
 * Precondition: kmerID_ trimmed to largest prefix encoded length
 * Output: kmerID trimmed to largest length bit
 */
void chemical_filter_single_pass(TKmerID & kmerID)
{
    // TKmerID code = kmerID & ~PREFIX_SELECTOR;
    auto [prefix, code] = split_kmerID(kmerID);
    assert(kmerID > 0 && code > 0);
    uint8_t AT = 0; // counter 'A'|'T'
    uint8_t CG = 0; // counter 'C'|'G'

    uint64_t enc_l = PRIMER_MAX_LEN - ffsll(prefix >> 54) + 1;
    // sum CG, AT content for longest kmer
    while (code != 1)
    {
        switch(3 & code)
        {
            case 1:
            case 2: ++CG; break;
            default: ++AT;
        }
        // TATA box test

        if ((code > (1ULL << 9)) && (((code & 0b11111111) == 0b11001100) || ((code & 0b11111111) == 0b00110011)))
        {
            if (enc_l <= PRIMER_MIN_LEN) // delete all length bits to discard this kmer
            {
                kmerID &= ~PREFIX_SELECTOR;
                return;
            }
            else // delete all length bits between current encoded length and larger ones
            {
                prefix &= ~((1ULL << (WORD_SIZE - enc_l + PRIMER_MIN_LEN)) - 1);
            }
        }
        if (!prefix)
            return;
        code >>= 2;
        --enc_l;
    }
    code = kmerID & ~PREFIX_SELECTOR;
    kmerID = prefix | code;
    // start with largest kmer
    uint64_t tailing_zeros = ffsll(prefix >> 54) - 1;
    uint64_t mask = 1ULL << (54 + tailing_zeros);
    for (uint8_t i = 0; i < 10 - tailing_zeros; ++i)
    {
        if (kmerID & mask)
        {
            // reset bit if Tm out of range
            auto Tm = (AT << 1) + (CG << 2);
            if (Tm < PRIMER_MIN_TM || Tm > PRIMER_MAX_TM)
            {
                kmerID ^= mask;
            }
            else
            {
                float CG_content = float(CG) / (float(__builtin_clzl(mask) + PRIMER_MIN_LEN));
                if (CG_content < CG_MIN_CONTENT || CG_content > CG_MAX_CONTENT)
                    kmerID ^= mask;
            }
        }
        switch(3 & code)
        {
            case 1:
            case 2: --CG; break;
            default: --AT;
        }
        code >>= 2;
        mask <<= 1;
    }

    // Filter di-nucleotide repeats and
    filter_repeats_runs(kmerID);
    // if (!(kmerID & PREFIX_SELECTOR) && passed)
    // {
    //     std::cout << "Did not pass filter_repeats_runs." << std::endl;
    //     passed = 0;
    // }
    // // Check for self-annealing
    // filter_annealing_connected(kmerID);
    // if (!(kmerID & PREFIX_SELECTOR) && passed)
    // {
    //     passed = 0;
    //     std::cout << "Did not pass filter_annealing_connected." << std::endl;
    // }
    // filter_annealing_disconnected(kmerID);
    // if (!(kmerID & PREFIX_SELECTOR) && passed)
    // {
    //     std::cout << "Did not pass filter_annealing_disconnected." << std::endl;
    // }

}

/*
 * Test for self-annealing.
 *
 */
/*void chemical_filter_square(TKmerID & kmerID)
{

}
*/

/* Helper function for computing the convolution of two sequences. For each overlap
 *position the Gibb's free energy is computed and the minimum returned;
 * TODO: see "Improved thermodynamic parameters and helix initiation factor to predict stability of DNA duplexes" Sugimoto et al. 1996
 */
extern inline float gibbs_free_energy(seqan::String<priset::dna> const & s, seqan::String<priset::dna> const & t)
{
    int8_t offset = 2;
    int8_t const n = seqan::length(s);
    int8_t const m = seqan::length(t);
    int8_t s1, s2;
    int8_t t1, t2;
    int8_t AT, CG;
    auto energy = [](int8_t const & CG, int8_t const &AT) { return -(CG + 2*AT);};
    //auto lambda = [](const std::string& s) { return std::stoi(s); };
    int8_t energy_min = 0;
    for (int8_t i = 0; i < n + m - 2 * offset + 1; ++i)
    {
        s1 = (n - offset - i < 0) ? 0 : n - offset - i; // std::max<int8_t>(n - offset - i, 0);
        t1 = std::max<int8_t>(i - n + offset, 0);
        t2 = std::min<int8_t>(i + offset, m);
        s2 = std::min<int8_t>(s1 + t2 - t1, n);
        // count complementary bps
        CG = 0;
        AT = 0;
        for (auto j = 0; j < s2-s1; ++j)
        {
            switch(char(s[s1+j]))
            {
                case 'A': AT += (t[t1+j] == 'T') ? 1 : 0; break;
                case 'T': AT += (t[t1+j] == 'A') ? 1 : 0; break;
                case 'C': CG += (t[t1+j] == 'G') ? 1 : 0; break;
                case 'G': CG += (t[t1+j] == 'C') ? 1 : 0; break;
                default: std::cout << "ERROR: primer contains unknown symbol '" << s[s1+j] << "'" << std::endl;
            }
            // update minimal energy
            if (energy(CG, AT) < energy_min)
                energy_min = energy(CG, AT);
        }
    }
    return energy_min;
}

} // namespace priset
