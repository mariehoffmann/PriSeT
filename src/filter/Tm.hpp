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

std::string kmerID2str(TKmerID kmerID);
extern inline void reset_length_leq(TKmerID & kmerID, uint8_t l);

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
extern inline void Tm_filter(TKmerID & kmer, uint8_t const Tm_min, uint8_t const Tm_max)
{
    uint64_t const l_sv = encoded_length(kmer);
    uint64_t const head_selector = 0b11ULL << 62;
    std::cout << "head_selector = " << std::bitset<64>(head_selector) << std::endl;
    std::cout << "kmer = " << kmerID2str(kmer) << std::endl;
    uint64_t l{0};
    auto [prefix, code] = split_kmerID(kmer);
    uint64_t code_sv{code};
    // std::cout << "length = " << l << endl;
    code <<= WORD_SIZE - l_sv + 2;
    // cout << "code = " << bitset<64>(code) << endl;
    uint8_t AT{0};
    uint8_t CG{0};
    while (l <= l_sv && prefix)
    {
        std::cout << "head = " << ((code & head_selector) >> 62) << std::endl;
        switch (code & head_selector)
        {
            case 0:
            case 3: AT += 2; break;
            case 1:
            case 2: CG += 4;
        }
        ++++l;
        std::cout << "l = " << int(l) << ", Tm = " << int(AT + CG) << ", AT = " << int(AT) << ", CG = " << int(CG) << std::endl;
        code <<= 2;
        if (l < (KAPPA_MIN << 1))
            continue;

        if (AT + CG > Tm_max) // invalidate all kmers larger or equal than length
        {
            std::cout << "reset all larger bits for l = " << int(l) << std::endl;
            kmer = prefix | code_sv;
            reset_length_leq(kmer, l);
            return;
        }
        uint64_t mask_l = (prefix & (1ULL << (WORD_SIZE - 1 - (l >> 1) + KAPPA_MIN)));
        if (mask_l > 0 && ((AT + CG) < Tm_min)) // invalidate only current kmer
        {
            prefix = reset_length(prefix, l);
            std::cout << "reset bit for l = " << int(l) << std::endl;
        }
    }
    kmer = prefix | code_sv;
}

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
    auto target_l = KAPPA_MIN;
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

// Difference in melting temperatures (degree Celsius) according to Wallace rule.
extern inline uint8_t dTm(TKmerID const kmerID1, TKmerID const mask1, TKmerID const kmerID2, TKmerID const mask2)
{
    // remove length mask
    TKmerID code1 = kmerID1 & ~PREFIX_SELECTOR;
    TKmerID code2 = kmerID2 & ~PREFIX_SELECTOR;

    uint8_t enc_l1 = encoded_length(code1); // encoded length
    uint8_t mask_l1 = (__builtin_clzl(mask1) + KAPPA_MIN) << 1; // selected length

    uint8_t enc_l2 = encoded_length(code2); // encoded length
    uint8_t mask_l2 = (__builtin_clzl(mask2) + KAPPA_MIN) << 1; // selected length

    // if (mask_l1 > enc_l1 || mask_l2 > enc_l2)
    //     std::cerr << "ERROR: largest encoded kmer length undershoots target length!\n";

    // trim kmer if exceeding encoded length
    code1 >>= enc_l1 - mask_l1;
    code2 >>= enc_l2 - mask_l2;

    int8_t AT{0};
    int8_t CG{0};
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

}
