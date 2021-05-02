// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================
//          Author: Marie Hoffmann <ozymandiaz147 AT gmail.com>
//          Manual: https://github.com/mariehoffmann/PriSeT

#pragma once

#include <iostream>
#include <fstream>
#include <numeric>
#include <vector>

#include "../common.hpp"
#include "annealing.hpp"
#include "CG.hpp"
#include "Tm.hpp"
#include "repeats.hpp"
#include "types/all.hpp"
#include "AT_tail.hpp"

namespace priset
{

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
void filter_Cs(TKmerID & kmerID, PrimerConfig const & primer_cfg)
{
    // melting temperature test
    Tm_filter(kmerID, primer_cfg.get_Tm_min(), primer_cfg.get_Tm_max(), primer_cfg.get_kappa_min(), primer_cfg.get_kappa_max());
    if (!(kmerID & PREFIX_SELECTOR))
        return;
    CG_filter(kmerID, primer_cfg.get_CG_min(), primer_cfg.get_CG_max(), primer_cfg.get_kappa_min(), primer_cfg.get_kappa_max());
    if (!(kmerID & PREFIX_SELECTOR))
        return;
    
    auto [prefix, code] = split_kmerID(kmerID);
    assert(kmerID > 0 && code > 0);
    uint8_t AT = 0; // counter 'A'|'T'
    uint8_t CG = 0; // counter 'C'|'G'

    // TATA box test
    uint64_t enc_l = KAPPA_MAX - ffsll(prefix >> 54) + 1;
    while (code != 1)
    {
        if ((code > (1ULL << 9)) && (((code & 0b11111111) == 0b11001100) || ((code & 0b11111111) == 0b00110011)))
        {
            if (enc_l <= KAPPA_MIN) // delete all length bits to discard this kmer
            {
                kmerID &= ~PREFIX_SELECTOR;
                return;
            }
            else // delete all length bits between current encoded length and larger ones
            {
                prefix &= ~((1ULL << (WORD_SIZE - enc_l + KAPPA_MIN)) - 1);
            }
        }
        if (!prefix)
            return;
        code >>= 2;
        --enc_l;
    }
    code = kmerID & ~PREFIX_SELECTOR;
    kmerID = prefix | code;

    filter_repeats_runs(kmerID);
    if (!(kmerID & PREFIX_SELECTOR))
        return;
    
    self_annealing_filter(kmerID);   
}

}
