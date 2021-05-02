// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================
//          Author: Marie Hoffmann <ozymandiaz147 AT gmail.com>
//          Manual: https://github.com/mariehoffmann/PriSeT

#pragma once

#include <numeric>
#include <vector>

#include "../types/IOConfig.hpp"
#include "../types/PrimerPair.hpp"
#include "../types/PrimerPairUnpacked.hpp"
#include "../types/PrimerConfig.hpp"
#include "utilities.hpp"

namespace priset
{

std::string kmerID2str(TKmerID kmerID);
extern inline void reset_length_leq(TKmerID & kmerID, uint8_t l);

// l1, l2 refer to kmer lengths including dimerizing 4-mer
extern inline bool annealing_helper(TKmerID & kmerID1, uint8_t l1, TKmerID & kmerID2, uint8_t l2)
{
    if (!kmerID2) // first 4-mer occurrence significant
        reset_length_leq(kmerID1, std::min(l1, l2));
    else
    {
        // case: l1 <= 16 and l2 > 16 or l1 > 16 and l2 > 16
        if (l2 > (KAPPA_MIN << 1) && l2 > l1)
        {
            reset_length_leq(kmerID2, l2);
        }
        // case: l1 > 16 and l2 <= 16 or l1 > 16 and l2 > 16
        if (l1 > (KAPPA_MIN << 1) && l1 > l2)
        {
            reset_length_leq(kmerID1, l1);
        }
        // case: l1 <= 16 and l2 <= 16
        if (l1 <= (KAPPA_MIN << 1) && l2 <= (KAPPA_MIN << 1))
        {
            reset_length_leq(kmerID1, l1);
            reset_length_leq(kmerID2, l2);
        }
    }
    return (!(kmerID1 & PREFIX_SELECTOR) || !(kmerID2 & PREFIX_SELECTOR)) ? 1 : 0;
}

/*
 * Self- and cross-dimerization test (if 2nd kmer is not zero). 4-mers ending before
 * position KAPPA_MIN result in deletion of length mask, otherwise only
 * affected lengths are erased. In the second part non-consecutive annealing is
 * tested - more than 8 annealing positions are considered to be critical.
 */
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
    for (overlap = 8; (overlap_mask = (1ULL << overlap) - 1) <= (1ULL << (l1 - 2)) - 1; ++++overlap)
    {
        // std::cout << "current overlap = " << int(overlap) << std::endl;
        if (ctr < 0)
            exit(0);
        --ctr;
        // prefix of code1 and suffix of code2
        prefix = (code1 >> (l1 - overlap));
        // handle case l1 > l2, trim leading bits from 1st code
        if (kmerID2 && overlap > l2 && l2 < l1)
            overlap_mask = (1 << l2) - 1;
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

    trim_to_true_length(kmerID1);
    if (kmerID2)
        trim_to_true_length(kmerID2);
    auto [prefix1, code1] = split_kmerID(kmerID1);
    auto [prefix2, code2] = (kmerID2) ? split_kmerID(kmerID2) : split_kmerID(kmerID1);
    uint64_t code2_rev = reverse(code2);
    uint8_t l1 =  WORD_SIZE - __builtin_clzll(code1) - 1;
    uint8_t l2 =  WORD_SIZE - __builtin_clzll(code2) - 1;

    uint8_t overlap = 20;
    uint64_t overlap_mask;
    uint64_t prefix, x, y;

    // Due to symmetry we need to test only the first half of possible overlapping positions
    for (overlap = 20; (overlap_mask = (1ULL << overlap) - 1) <= (1ULL << (std::min(l1, l2) - 2)) - 1; ++++overlap)
    {
        // prefix of code1 and suffix of code2
        prefix = (code1 >> (l1 - overlap));
        // creates 11 patterns with even offsets for complementary bases
        x = (prefix ^ code2) & overlap_mask;
        y = (prefix ^ code2_rev) & overlap_mask;
        // count ones every two bits and store in 2 bits, 011011 -> 010110
        x = x - ((x >> 1) & 0x1555555555555);
        // filter every second set bit
        x &= 0xAAAAAAAAAAAA & overlap_mask; // =0b10_{50}
        auto annealings_x = __builtin_popcountll(x);
        if (annealings_x >= 8)
            annealing_helper(kmerID1, overlap, kmerID2, l2);
        y = y - ((y >> 1) & 0x1555555555555);
        // filter every second set bit
        y &= 0xAAAAAAAAAAAA; // =0b10_{42}
        auto annealings_y = __builtin_popcountll(y);
        if (annealings_y >= 8)
            annealing_helper(kmerID1, overlap, kmerID2, l2);
    }

    if (kmerID1 & PREFIX_SELECTOR)
    {
        overlap = (l1 <= l2) ? l1 : l2;
        while (overlap >= 20)
        {
            overlap_mask = (1ULL << overlap) - 1;
            if (kmerID2) // cross-annealing with code2 forward
            {
                x = (code1 ^ code2) & overlap_mask;
                x = x - ((x >> 1) & 0x1555555555555);
                x &= 0xAAAAAAAAAAAA & overlap_mask; // =0b10_{50}
                auto annealings_x = __builtin_popcountll(x);
                if (annealings_x >= 8)
                {
                    kmerID1 &= ~PREFIX_SELECTOR;
                    kmerID2 &= ~PREFIX_SELECTOR;
                    break;
                }
            }
            y = (code1 ^ code2_rev) & overlap_mask;
            y = y - ((y >> 1) & 0x1555555555555);
            y &= 0xAAAAAAAAAAAA & overlap_mask; // =0b10_{42}
            auto annealings_y = __builtin_popcountll(y);
            if (annealings_y >= 8)
            {
                kmerID1 &= ~PREFIX_SELECTOR;
                kmerID2 &= ~PREFIX_SELECTOR;
                break;
            }
            code2 >>= 2;
            overlap = std::min(WORD_SIZE - __builtin_clzll(code2) - 1, int(l1));
        }
    }
}

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
            if (l1 <= KAPPA_MIN && l2 <= KAPPA_MIN)
            {
                kmerID1 &= ~PREFIX_SELECTOR;
                kmerID2 &= ~PREFIX_SELECTOR;
                return;
            }
            else
            {
                // kmerID1 remains combinable with prefix of kmerID2 or l2 prefix needs to be trimmed
                if (l1 <= KAPPA_MIN || l2 > KAPPA_MIN)
                    kmerID2 = (kmerID2 & ~PREFIX_SELECTOR) + (prefix2 & ~((ONE_LSHIFT_63 >> (l2 - KAPPA_MIN - 1)) - 1));
                if (l2 <= KAPPA_MIN || l1 > KAPPA_MIN) // kmerID2 remains combinable with prefix of kmerID1
                    kmerID1 = (kmerID1 & ~PREFIX_SELECTOR) + (prefix1 & ~((ONE_LSHIFT_63 >> (l1 - KAPPA_MIN - 1)) - 1));
            }
        }
        code2 >>= 2;
    }
}

/*
* Helper to check for self-annealing patterns (connected and disconnected).
*/
void extern inline self_annealing_filter(TKmerID & kmerID)
{
    filter_annealing_connected(kmerID);
    if (kmerID & PREFIX_SELECTOR)
        filter_annealing_disconnected(kmerID);
}

void extern inline cross_annealing_filter(TKmerID & kmerID1, TKmerID & kmerID2)
{
    filter_annealing_connected(kmerID1, kmerID2);
    if ((kmerID1 & PREFIX_SELECTOR) && (kmerID2 & PREFIX_SELECTOR))
        filter_annealing_disconnected(kmerID1, kmerID2);
}

}
