#pragma once

#include <numeric>
#include <vector>

#include "types/IOConfig.hpp"
#include "types/PrimerPair.hpp"
#include "types/PrimerPairUnpacked.hpp"
#include "types/PrimerConfig.hpp"
#include "utilities.hpp"

namespace priset
{

/*
 * Detect and filter mononucleotide runs and dinucleotide repeats of length four,
 * and five, respectively.
 */
extern inline void filter_repeats_runs(TKmerID & kmerID)
{
    auto [prefix, code] = split_kmerID(kmerID);
    uint64_t const tail_selector_10 = (1 << 10) - 1;
    uint64_t const tail_selector_20 = (1 << 20) - 1;
    kmerID = code; // save trimmed kmer-code part
    uint64_t const k = (63 - __builtin_clzl(code)) >> 1;
    uint64_t offset;
    for (uint64_t i = 0; i < k - 4; ++i) // up-to four, i.e. 8 bits consecutive characters permitted
    {
        uint64_t tail_10 = tail_selector_10 & code;
        // XOR tail with A_4, C_4, G_4, T_4
        if (    (tail_10 == 0b00000000) || (tail_10 == 0b01010101) ||
                (tail_10 == 0b10101010) || (tail_10 == 0b11111111))
        {
            // note: (x >> 64) won't be executed
            offset = 64 - (std::max((uint64_t)KAPPA_MIN, k - i) - KAPPA_MIN);
            // delete all k bits ≥ len_selector
            prefix = (offset == 64) ? 0 : (prefix >> offset) << offset;
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
                offset = 64 - (std::max((uint64_t)KAPPA_MIN, k - i) - KAPPA_MIN);
                // delete all k bits ≥ len_selector
                prefix = (offset == 64) ? 0 : (prefix >> offset) << offset;
            }
        }
        if (!prefix)
            break;
        code >>= 2;
    }
    // add updated prefix
    kmerID |= prefix;
}

}
