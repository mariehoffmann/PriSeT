// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================
//          Author: Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
//          Manual: https://github.com/mariehoffmann/PriSeT

// Global Functions for Primer Constraint Checking.

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <functional>
#include <vector>

#include <seqan/basic.h>

#include "primer_cfg_type.hpp"
#include "types.hpp"
#include "utilities.hpp"

namespace priset  //::chemistry TODO: introduce chemistry namespace
{

//!\brief Wallace rule to compute the melting temperature of a primer sequence given as 64 bit code.
float primer_melt_wallace(TKmerID code)
{
    uint8_t ctr_AT = 0, ctr_CG = 0;
    std::array<char, 4> decodes = {'A', 'C', 'G', 'T'};
    while (code != 1)
    {
        switch(decodes[3 & code])
        {
            case 'A':
            case 'T': ++ctr_AT; break;
            case 'C':
            case 'G': ++ctr_CG;
        }
        code >>= 2;
    }
    return (ctr_AT << 1) + (ctr_CG << 2);
}

//!\brief Salt-adjusted method to compute melting temperature of primer sequence.
// input primer:string sequence, Na:float molar Natrium ion concentration
inline float primer_melt_salt(TKmerID code, float const Na)
{
    uint8_t ctr_CG = 0;
    uint8_t seq_len = 0;
    std::array<char, 4> decodes = {'A', 'C', 'G', 'T'};
    while (code != 1)
    {
        switch(decodes[3 & code])
        {
            case 'C':
            case 'G': ++ctr_CG;
        }
        code >>= 2;
        ++seq_len;
    }
    return 100.5 + 41.0 * ctr_CG / seq_len - 820.0 / seq_len + 16.6 * std::log10(Na);
}

//!\brief Compute melting temperature of primer sequence with method set in primer configuration.
float get_Tm(primer_cfg_type const & primer_cfg, TKmerID kmer_ID) noexcept
{
    switch(primer_cfg.get_primer_melt_method())
    {
        case TMeltMethod::WALLACE: return primer_melt_wallace(kmer_ID);
        default: return primer_melt_salt(kmer_ID, primer_cfg.get_Na());
    }
}

//!\brief Check if CG content is in range set by the primer configurator.
// Returns false if constraint is violated.
bool filter_CG(primer_cfg_type const & primer_cfg, TKmerID code)
{
    uint8_t ctr_CG = 0;
    uint8_t seq_len = 0;
    std::array<char, 4> decodes = {'A', 'C', 'G', 'T'};
    while (code != 1)
    {
        switch(decodes[3 & code])
        {
            case 'C':
            case 'G': ++ctr_CG; break;
        }
        code >>= 2;
        ++seq_len;
    }
    // relative CG content
    float CG = float(ctr_CG) / float(seq_len);
    return (CG >= primer_cfg.get_min_CG_content() && CG <= primer_cfg.get_max_CG_content());
}

//  Check if not more than 3 out of the 5 last bases at the 3' end are CG.
//  DNA sense/'+': 5' to 3', antisense/'-': 3' to 5'
// Returns false if constraint is violated.
inline bool filter_CG_clamp(/*primer_cfg_type const & primer_cfg, */seqan::String<priset::dna> const & sequence, char const sense)
{
    assert(seqan::length(sequence) < (1 << 8));
    assert(sense == '+' || sense == '-');
    uint8_t CG_cnt = 0;
    uint8_t offset = (sense == '+') ? 0 : seqan::length(sequence) - 6;
    for (uint8_t i = 0; i < 5; ++i)
        CG_cnt += (sequence[i + offset] == 'C' || sequence[i + offset] == 'G') ? 1 : 0;
    return CG_cnt <= 3;
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

/*
 * Filter di-nucleotide repeats, like ATATATAT, and runs, i.e. series of same nucleotides.
 * For both the maximum is 4 consecutive di-nucleotides, and 4bp, respectively.
 * Needs to be called with non-ambigous kmerID, i.e. an ID encoding a length-variable kmer sequence.
 */
bool filter_repeats_runs(TKmerID kmer_ID)
{
    TSeq seq = dna_decoder(kmer_ID);
    if (seqan::length(seq) > 4)
    {
        uint8_t repeat_even = 1;
        uint8_t repeat_odd = 1;
        //uint8_t & repeat;
        TSeq ifx_even = seqan::infixWithLength(seq, 0, 2); // even start positions
        TSeq ifx_odd = seqan::infixWithLength(seq, 1, 2); // odd start positions
        //TSeq & ifx;
        TSeq aux;
        for (uint8_t i = 3; i < seqan::length(seq); ++i)
        {
            aux = seqan::infixWithLength(seq, i - 1, 2);
            TSeq & ifx = (i % 2) ? ifx_even : ifx_odd;
            uint8_t & repeat = (i % 2) ? repeat_even : repeat_odd;
            if (aux[0] == ifx[0] && aux[1] == ifx[1])
            {
                ++repeat;
                if (repeat >= 4)
                    return false;
            }
            else
            {
                repeat = 1;
                ifx = aux;
            }
        }

        seqan::Finder<TSeq> finder(seq);
        for (char c : std::vector<char>{'A', 'C', 'G', 'T'})
        {
            seqan::Pattern<TSeq, Horspool> pattern(std::string(5, c));
            if (seqan::find(finder, pattern))
                return false;
        }
    }
    return true;
}

/* Helper function for computing the convolution of two sequences. For each overlap
 *position the Gibb's free energy is computed and the minimum returned;
 * TODO: see "Improved thermodynamic parameters and helix initiation factor to predict stability of DNA duplexes" Sugimoto et al. 1996
 */
float gibbs_free_energy(seqan::String<priset::dna> const & s, seqan::String<priset::dna> const & t)
{
    int8_t offset = 2;
    int8_t const n = seqan::length(s);
    int8_t const m = seqan::length(t);
    int8_t s1, s2;
    int8_t t1, t2;
    int8_t ctr_AT, ctr_CG;
    auto energy = [](int8_t const & ctr_CG, int8_t const &ctr_AT) { return -(ctr_CG + 2*ctr_AT);};
    //auto lambda = [](const std::string& s) { return std::stoi(s); };
    int8_t energy_min = 0;
    for (int8_t i = 0; i < n + m - 2 * offset + 1; ++i)
    {
        s1 = (n - offset - i < 0) ? 0 : n - offset - i; // std::max<int8_t>(n - offset - i, 0);
        t1 = std::max<int8_t>(i - n + offset, 0);
        t2 = std::min<int8_t>(i + offset, m);
        s2 = std::min<int8_t>(s1 + t2 - t1, n);
        // count complementary bps
        ctr_CG = 0;
        ctr_AT = 0;
        for (auto j = 0; j < s2-s1; ++j)
        {
            switch(char(s[s1+j]))
            {
                case 'A': ctr_AT += (t[t1+j] == 'T') ? 1 : 0; break;
                case 'T': ctr_AT += (t[t1+j] == 'A') ? 1 : 0; break;
                case 'C': ctr_CG += (t[t1+j] == 'G') ? 1 : 0; break;
                case 'G': ctr_CG += (t[t1+j] == 'C') ? 1 : 0; break;
                default: std::cout << "ERROR: primer contains unknown symbol '" << s[s1+j] << "'" << std::endl;
            }
            // update minimal energy
            if (energy(ctr_CG, ctr_AT) < energy_min)
                energy_min = energy(ctr_CG, ctr_AT);
        }
    }
    return energy_min;
}

/* !\brief Check for self-dimerization, i.e. bonding energy by same sense bonding.
 * Returns true if ΔG ≥ -5kcal/mol
 */
inline bool filter_cross_dimerization(TKmerID kmer_ID1, TKmerID kmer_ID2)
{
//    std::cout << "enter filter_cross_dimerization with kmer_ID1 = " << kmer_ID1 << " and kmer_ID2 = " << kmer_ID2 << std::endl;
    std::vector<TSeq> decodes;
    dna_decoder(kmer_ID1, decodes);
    dna_decoder(kmer_ID2, decodes);
    auto seq1 = decodes.at(0);
    auto seq2 = decodes.at(1);
    float dG = gibbs_free_energy(seq1, seq2);
    //std::cout << "minimal free energy for self-dimerization of s,t is " << dG << std::endl;
    return (dG < -10) ? false : true;
}

/* !\brief Check for self-dimerization, i.e. bonding energy by same sense bonding.
 * Returns true if ΔG ≥ -5kcal/mol
 */
inline bool filter_self_dimerization(TKmerID kmer_ID)
{
    return filter_cross_dimerization(kmer_ID, kmer_ID);
}

/*

# check for 2ndary structure hairpin, may only be present at 3' end with a delta(G) = -2 kcal/mol,
# or internally with a delta(G) of -3 kcal/mol
# TODO: upper limit for loop length is disrespected currently
def filter_hairpin(seq, cfg):
    n, min_loop_len = len(seq), int(cfg.var['hairpin_loop_len'][0])
    palindrome_len_rng = range(3, len(seq)/2 - min_loop_len + 1)
    seq_ci = complement(seq)[::-1] # inverted, complemented sequence
    for i in range(len(seq) - 2*palindrome_len_rng[0] - min_loop_len):
        for m in palindrome_len_rng:
            for i_inv in range(n - 2*m - min_loop_len):
                if seq[i:i+m] == seq_ci[i_inv:i_inv+m]:
                    #print seq[i:i+m], ' == ', seq_ci[i_inv:i_inv+m]
                    return False
    return True

'''
    Same sense interaction: sequence is partially homologous to itself.
'''
def filter_selfdimer(seq, cfg):
    seq_rev = seq[::-1]
    return filter_crossdimer(seq, seq[::-1], cfg)

'''
    Pairwise interaction: sequence s is partially homologous to sequence t.
    If primer binding is too strong in terms of delta G, less primers are available for DNA binding.
    Todo: Use discrete FFT convolution to compute all free Gibb's energy values for all overlaps in O(nlogn).
'''
def filter_crossdimer(s, t, cfg):
    n, m = len(s), len(t)
    cnv = [0 for _ in range(len(s)+len(t)-1)]
    for n in range(len(s) + len(t) - 1):
        cnv[n] = sum([delta_G(s[m], t[n-m]) for m in range(len(s))])
    if min(cnv) < cfg.var['delta_G_cross']:
        return False
    return True, min(cnv)
*/

} // namespace chemistry
