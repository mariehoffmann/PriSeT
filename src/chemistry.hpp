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

std::string dna_decoder(uint64_t code, uint64_t const mask);

// Difference in melting temperatures (degree Celsius) according to Wallace rule.
extern inline float dTm(TKmerID const kmerID1_, TKmerID const mask1, TKmerID const kmerID2_, TKmerID const mask2)
{
    if (!(kmerID1_ & mask1) || !(kmerID2_ & mask2))
    {
        std::cout << "Error: target mask bit not set\n";
        exit(0);
    }
    // remove length mask
    TKmerID kmerID1 = (kmerID1_ << LEN_MASK_SIZE) >> LEN_MASK_SIZE;
    TKmerID kmerID2 = (kmerID2_ << LEN_MASK_SIZE) >> LEN_MASK_SIZE;

    auto enc_l1 = (WORD_SIZE - 1 - __builtin_clzl(kmerID1)) >> 1; // encoded length
    auto mask_l1 = __builtin_clzl(mask1) + PRIMER_MIN_LEN;      // selected length
    kmerID1 >>= (enc_l1 - mask_l1) << 1;   // string correction

    auto enc_l2 = (WORD_SIZE - 1 - __builtin_clzl(kmerID2)) >> 1; // encoded length
    auto mask_l2 = __builtin_clzl(mask2) + PRIMER_MIN_LEN;      // selected length
    kmerID2 >>= (enc_l2 - mask_l2) << 1;   // string correction

    //std::cout << "corrected kmerID1 = " << dna_decoder(kmerID1) << std::endl;
    //std::cout << "corrected kmerID2 = " << dna_decoder(kmerID2) << std::endl;
    int8_t ctr_AT = 0;
    int8_t ctr_CG = 0;
    if (!(__builtin_clzl(kmerID1) % 2))
    {
        std::cout << "ERROR: expected odd number of leading zeros\n";
        exit(0);
    }
    if (!(__builtin_clzl(kmerID2) % 2))
    {
        std::cout << "ERROR: expected odd number of leading zeros\n";
        exit(0);
    }
    while (kmerID1 != 1)
    {
        if (!(kmerID1 & 3) || (kmerID1 & 3) == 3)  // 'A' (00) or 'T' (11)
            ++ctr_AT;
        else
            ++ctr_CG;
        kmerID1 >>= 2;
    }
    while (kmerID2 != 1)
    {
        if (!(kmerID2 & 3) || (kmerID2 & 3) == 3)  // 'A' (00) or 'T' (11)
            --ctr_AT;
        else
            --ctr_CG;
        kmerID2 >>= 2;
    }
    return (std::abs(ctr_AT) << 1) + (abs(ctr_CG) << 2);
}

//TAACGTAACGTACGTACGT, TAGCTAACTACATAGCTACGA
// AT = 11              13
// CG =  8              8

//!\brief Wallace rule to compute the melting temperature of a primer sequence given as 64 bit code.
extern inline float primer_melt_wallace(TKmerID code)
{
    uint8_t ctr_AT = 0, ctr_CG = 0;
    std::array<char, 4> decodes = {'A', 'C', 'G', 'T'};
    code &= ~((1ULL << 52ULL) - 1ULL); // delete length mask
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
extern inline float primer_melt_salt(TKmerID code, float const Na)
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
extern inline float get_Tm(primer_cfg_type const & primer_cfg, TKmerID kmer_ID) noexcept
{
    switch(primer_cfg.get_primer_melt_method())
    {
        case TMeltMethod::WALLACE: return primer_melt_wallace(kmer_ID);
        default: return primer_melt_salt(kmer_ID, primer_cfg.get_Na());
    }
}

//!\brief Check if CG content is in range set by the primer configurator.
// Returns false if constraint is violated.
extern inline bool filter_CG(primer_cfg_type const & primer_cfg, TKmerID code)
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
extern inline bool filter_CG_clamp(/*primer_cfg_type const & primer_cfg, */seqan::String<priset::dna> const & sequence, char const sense)
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
extern inline bool filter_repeats_runs(TKmerID kmerID)
{
    for (uint64_t mask = ONE_LSHIFT_63; mask >= (1ULL << 53); mask >>= 1)
    {
        if (mask & kmerID)
        {
            TSeq seq = dna_decoder(kmerID, mask);
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
                        if (repeat > 4)
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
        }
    }
    return true;
}

uint64_t code_prefix(uint64_t const code_, uint64_t mask);

template<typename uint_type>
std::string bits2str(uint_type i);

std::string code2str(TKmerID kmerID);

extern inline void filter_repeats_runs2(TKmerID & kmerID_)
{
    //std::cout << "Enter filter_repeats_runs2\n";
    uint64_t mask = kmerID_ & MASK_SELECTOR;
    uint64_t const tail_selector_10 = (1 << 10) - 1;
    uint64_t const tail_selector_20 = (1 << 20) - 1;
    TKmerID suffix = code_prefix(kmerID_, 0); // trim to true length
    std::cout << "true length trimmed: " << bits2str(suffix) <<std::endl;
    kmerID_ = suffix; // save trimmed kmer-code part
    //std::cout << "head-less sequence: " << code2str(kmerID) << " and head = " << bits2str(mask>>54) << std::endl;
    uint64_t const k = (63 - __builtin_clzl(suffix)) >> 1;
    for (uint64_t i = 0; i < k - 4; ++i) // up-to four, i.e. 8 bits consecutive characters permitted
    {
        uint64_t tail_10 = tail_selector_10 & suffix;
        // XOR tail with A_5, C_5, G_5, T_5
        if ((tail_10 == 0b0000000000) || (tail_10 == 0b0101010101) || (tail_10 == 0b1010101010) || (tail_10 == 0b1111111111))
        {
            // note that (x >> 64) won't be executed
            uint64_t offset = 64 - (std::max(PRIMER_MIN_LEN, k - i) - PRIMER_MIN_LEN);
            // delete all k bits ≥ len_selector
            mask = (offset == 64) ? 0 : (mask >> offset) << offset;
        }
        // at least 10 characters left to test for di-nucleotide repeats of length 5
        if (k - i > 9)
        {
            uint64_t tail_20 = suffix & tail_selector_20;
            // AT_5, TA_5, AC_5, CA_5, AG_5, GA_5, CG_5, GC_5, CT_5, TC_5, GT_5, TG_5
            if ( (tail_20 == 0b00110011001100110011) || (tail_20 == 0b11001100110011001100) ||
                 (tail_20 == 0b00010001000100010001) || (tail_20 == 0b01000100010001000100) ||
                 (tail_20 == 0b00100010001000100010) || (tail_20 == 0b10001000100010001000) ||
                 (tail_20 == 0b01100110011001100110) || (tail_20 == 0b10011001100110011001) ||
                 (tail_20 == 0b01110111011101110111) || (tail_20 == 0b11011101110111011101) ||
                 (tail_20 == 0b10111011101110111011) || (tail_20 == 0b11101110111011101110))
            {
                //std::cout << "DEBUG: match for di-nucl run in tail: " << dna_decoder(suffix, 0) << std::endl;
                // delete all k bits ≥ len_selector
                uint64_t offset = 64 - (std::max(PRIMER_MIN_LEN, k - i) - PRIMER_MIN_LEN);
                // delete all k bits ≥ len_selector
                mask = (offset == 64) ? 0 : (mask >> offset) << offset;
                //std::cout << "DEBUG: new bit mask = " << bits2str<uint64_t>(mask >> 54) << std::endl;
            }
        }
        if (!mask)
            break;
        suffix >>= 2;
    }
    // add updated prefix
    kmerID_ += mask;
}

/*
 * Filter kmers based on their chemical properties regardless of their pairing.
 * This is a metafunction performing several single pass checks: melting tempaerature,
 * CG content, di-nucleotide repeats and consecutive runs.
 * Length bits in prefix of kmerID are reset in case of not passing.
 * This function only modifies the prefix. Precondition is that kmerID_ is trimmed,
 * i.e. largest prefix bit corresponds to encoded string length. The last filter trims
 * again to the maximal encoded length.
 */
void chemical_filter_single_pass(primer_cfg_type const & primer_cfg, TKmerID & kmerID_)
{
    assert(kmerID > 0);
    TKmerID code = (kmerID_ << LEN_MASK_SIZE) >> LEN_MASK_SIZE;
    float Tm = get_Tm(primer_cfg, kmerID);
    uint8_t AT_ctr = 0;
    uint8_t CG_ctr = 0;
    uint64_t length_bit_selector = ONE_LSHIFT_63;
    std::array<char, 4> decodes = {'A', 'C', 'G', 'T'};
    while (kmerID != 1)
    {
        switch(decodes[3 & code])
        {
            case 'C':
            case 'G': ++CG_ctr; break;
            default: ++AT_ctr;
        }
        if (kmerID_ & length_bit_selector)
        {
            // reset bit if Tm out of range
            if (((AT_ctr << 1) + (CG_ctr << 2) < PRIMER_MIN_TM) || ((AT_ctr << 1) + (CG_ctr << 2) > PRIMER_MAX_TM))
                kmerID_ ^= length_bit_selector;
            else
            {
                // reset bit if CG content out of range
                float CG = float(ctr_CG) / float(__builtin_clzl(length_bit_selector) + PRIMER_MIN_LEN);
                if (CG < CG_MIN_CONTENT || CG > CG_MAX_CONTENT)
                    kmerID_ ^= length_bit_selector;
            }
        }
        code >>= 2;
        length_bit_selector >>= 1;
    }
    // Filter di-nucleotide repeats and
    filter_repeats_runs2(kmerID_);
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
extern inline bool filter_cross_dimerization(TKmerID kmerID1, TKmerID kmerID2)
{
//    std::cout << "enter filter_cross_dimerization with kmerID1 = " << kmerID1 << " and kmerID2 = " << kmerID2 << std::endl;
    // Approximate: check for smallest encoded kmer
    // TODO: do check exact and avoid String conversion
    TSeq seq1 = dna_decoder(kmerID1, __builtin_clzl(kmerID1) + PRIMER_MIN_LEN);
    TSeq seq2 = dna_decoder(kmerID2, __builtin_clzl(kmerID2) + PRIMER_MIN_LEN);

    float dG = gibbs_free_energy(seq1, seq2);
    //std::cout << "minimal free energy for self-dimerization of s,t is " << dG << std::endl;
    return (dG < -10) ? false : true;
}

/* !\brief Check for self-dimerization, i.e. bonding energy by same sense bonding.
 * Returns true if ΔG ≥ -5kcal/mol
 */
extern inline bool filter_self_dimerization(TKmerID kmer_ID)
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
