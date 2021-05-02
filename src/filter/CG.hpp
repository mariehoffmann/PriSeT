// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================
//          Author: Marie Hoffmann <ozymandiaz147 AT gmail.com>
//          Manual: https://github.com/mariehoffmann/PriSeT

#pragma once

#include <numeric>
#include <vector>

#include "utilities.hpp"


namespace priset
{

std::string kmerID2str(TKmerID kmerID);
extern inline void reset_length_leq(TKmerID & kmerID, uint8_t l);

/*
 * Compute melting temperature of a single encoded oligomer. The oligomer is trimmed
 * to the length indicated by bit mask, otherwise (mask = 0) to the smallest encoded
 * length.
 * \Details
 * We use bit-parallelism to .
 */
extern inline float CG(TKmerID const kmerID, uint64_t const mask)
{
    auto [prefix, code] = split_kmerID(kmerID);
    uint8_t enc_l = encoded_length(kmerID);
    uint8_t target_l = KAPPA_MIN;
    target_l += (prefix & !mask) ? __builtin_clzll(prefix) : __builtin_clzll(mask);
    code >>= (enc_l - (target_l << 1));
    uint64_t p = 0x5555555555555ULL; // = (01)_52
    uint64_t q = 0xAAAAAAAAAAAAAULL; // = (10)_52 
    uint64_t x = ((code & p) << 1) ^ (code & q);
    return __builtin_popcountll(x) - 1;
}

extern inline float CG_percent(TKmerID const kmerID, uint64_t const mask)
{
    uint8_t target_l = KAPPA_MIN;
    auto [prefix, code] = split_kmerID(kmerID);
    target_l += (prefix & !mask) ? __builtin_clzll(prefix) : __builtin_clzll(mask);
    return float(CG(kmerID, mask))/float(target_l);
}

/*
 * Filter for encoded kmers to satisfy CG content.
 * K-mers that do not fail the CG constraint are deleted by reseting the
 * corresponding length bit.
 *
 * \Details
 * Compute the melting temperature applying the simple Wallace rule:
 *                            Tm = 2AT + 4CG
 * The difference in temperature estimates between the Wallace and the nearest
 * neighbour method (considered as more precise) rarely exceeds 4.5 K. This
 * maximum estimatation error needs to be considered when filtering and analyzing
 * potential primer pairs.
 */
extern inline void CG_filter(TKmerID & kmerID, float const CG_min, float const CG_max, uint8_t const kappa_min, uint8_t const kappa_max)
{
    uint64_t mask = (1ULL << 63) >> (kappa_min - KAPPA_MIN);
    float CG_content;
    for (uint8_t k = kappa_min; k <= kappa_max; ++k)
    {
        if (mask & kmerID)
        {
            CG_content = CG_percent(kmerID, mask);
            
            if (CG_content < CG_min || CG_content > CG_max)
                kmerID ^= mask;
        }
        mask >>= 1;
    }
}

}
