// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================
//          Author: Marie Hoffmann <ozymandiaz147 AT gmail.com>
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
#include "filter/CG.hpp"
#include "filter/Tm.hpp"
#include "types/all.hpp"
#include "utilities.hpp"

namespace priset  //::chemistry TODO: introduce chemistry namespace
{

std::string dna_decoder(uint64_t code, uint64_t const mask);
// std::pair<uint64_t, uint64_t> split(TKmerID kmerID);
extern inline uint64_t complement(uint64_t const code);
extern inline uint64_t reverse(uint64_t const code);
extern inline uint64_t reverse_complement(uint64_t const code);
extern inline void Tm_filter(TKmerID & kmerID, uint8_t const Tm_min, uint8_t const Tm_max, uint8_t const kappa_min, uint8_t const kappa_max);
extern inline void trim_to_true_length(TKmerID &);


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
/*bool filter_hairpin(PrimerConfig const & primer_cfg, seqan::String<priset::dna> const & sequence)
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
            uint8_t target_l = (KAPPA_MIN + __builtin_clzl(mask)) << 1; // in bits
            code >>= enc_l - target_l; // delete prefix and trim to target length
        }
    }
    return  ((((code & 3) | ((code & 3) + 1)) == 3) +
            ((((code >> 2) & 3) | (((code >> 2) & 3) + 1)) == 3) +
            ((((code >> 4) & 3) | (((code >> 4) & 3) + 1)) == 3) +
            ((((code >> 6) & 3) | (((code >> 6) & 3) + 1)) == 3) +
            ((((code >> 8) & 3) | (((code >> 8) & 3) + 1)) == 3) <= 3) ? true : false;
}


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
