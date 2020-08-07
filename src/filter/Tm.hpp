#pragma once

#include <numeric>
#include <vector>

#include "filter/CG.hpp"
#include "types/IOConfig.hpp"
#include "types/PrimerPair.hpp"
#include "types/PrimerPairUnpacked.hpp"
#include "types/PrimerConfig.hpp"
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
    uint8_t target_l = KAPPA_MIN;
    target_l += (prefix & !mask) ? __builtin_clzll(prefix) : __builtin_clzll(mask);
    uint8_t CG_content = CG(kmerID, mask);
    return ((target_l - CG_content) << 1) + (CG_content << 2);

    // code >>= (enc_l - (target_l << 1));
    // uint64_t p = 0x5555555555555ULL;
    // uint64_t q = 0xAAAAAAAAAAAAAULL;
    // uint64_t x = ((code & p) << 1) ^ (code & q);
    // uint8_t CG = __builtin_popcountll(x) - 1;
    // return ((target_l - CG) << 1) + (CG << 2);
}

/*
 * Filter encoded kmers for satisfying melting temperature (Tm) ranges.
 * K-mers that do not fail the Tm constraint are deleted by reseting the
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
extern inline void Tm_filter(TKmerID & kmerID, uint8_t const Tm_min, uint8_t const Tm_max, uint8_t const kappa_min, uint8_t const kappa_max)
{
    uint64_t mask = 1ULL << 63;
    uint8_t Tm_wallace;
    for (uint8_t k = kappa_min; k <= kappa_max; ++k)
    {
        if (mask & kmerID)
        {
            Tm_wallace = Tm(kmerID, mask);
            if (Tm_wallace < Tm_min)
                kmerID ^= mask;
            else if (Tm_wallace > Tm_max)
            {
                reset_length_leq(kmerID, encoded_length_mask(mask));
                return;
            }
        }
        mask >>= 1;
    }
}

// Difference in melting temperatures (degree Celsius) according to Wallace rule.
extern inline uint8_t dTm(TKmerID const kmerID1, TKmerID const mask1, TKmerID const kmerID2, TKmerID const mask2)
{
    return std::abs(Tm(kmerID1, mask1) - Tm(kmerID2, mask2));
}

}
