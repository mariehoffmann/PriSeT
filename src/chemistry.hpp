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

#include "primer_cfg_type.hpp"
#include "types.hpp"
#include "utilities.hpp"

namespace priset  //::chemistry TODO: introduce chemistry namespace
{

std::string dna_decoder(uint64_t code, uint64_t const mask);
std::pair<uint64_t, uint64_t> split(TKmerID kmerID);


// Difference in melting temperatures (degree Celsius) according to Wallace rule.
extern inline float dTm(TKmerID const kmerID1, TKmerID const mask1, TKmerID const kmerID2, TKmerID const mask2)
{
    if (!(kmerID1 & mask1) || !(kmerID2 & mask2))
        std::cerr << "Error: target mask bit not set\n";

    // remove length mask
    TKmerID code1 = kmerID1 & ~PREFIX_SELECTOR ;
    TKmerID code2 = kmerID2 & ~PREFIX_SELECTOR;

    uint8_t enc_l1 = (WORD_SIZE - 1 - __builtin_clzl(code1)) >> 1; // encoded length
    uint8_t mask_l1 = __builtin_clzl(mask1) + PRIMER_MIN_LEN;      // selected length

    uint8_t enc_l2 = (WORD_SIZE - 1 - __builtin_clzl(code2)) >> 1; // encoded length
    uint8_t mask_l2 = __builtin_clzl(mask2) + PRIMER_MIN_LEN;      // selected length

    if (mask_l1 > enc_l1 || mask_l2 > enc_l2)
        std::cerr << "ERROR: largest encoded kmer length undershoots target length!\n";

    // trim kmer if overlong
    code1 >>= (enc_l1 - mask_l1) << 1;
    code2 >>= (enc_l2 - mask_l2) << 1;

    //std::cout << "dTm: enc_l1 = " << enc_l1 << ", mask_l1 = " << mask_l1 << std::endl;
    //std::cout << "dTm: enc_l2 = " << enc_l2 << ", mask_l2 = " << mask_l2 << std::endl;

    //std::cout << "corrected code = " << dna_decoder(code1) << std::endl;
    //std::cout << "corrected code2 = " << dna_decoder(code2) << std::endl;
    int8_t AT = 0;
    int8_t CG = 0;

    if (!(__builtin_clzl(code1) % 2))
        std::cerr << "ERROR: expected odd number of leading zeros\n";

    if (!(__builtin_clzl(code2) % 2))
        std::cerr << "ERROR: expected odd number of leading zeros\n";

    while (code1 != 1)
    {
        // 'A' (00) or 'T' (11)
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

// Mumur formular Tm = 2AT + 4CG
// If no mask is given, the code is trimmed to the smallest encoded length.
extern inline uint8_t Tm(TKmerID const kmerID, uint64_t const mask)
{
    auto [prefix, code] = split(kmerID);
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
    return (std::abs(AT) << 1) + (abs(CG) << 2);
}

// CG content computation for output
extern inline float CG(TKmerID const kmerID, uint64_t const mask)
{
    auto [prefix, code] = split(kmerID);
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
    auto [prefix, code] = split(kmerID);
    if (!code)
        throw std::invalid_argument("Expected kmerID non zero!");
    if (!prefix)
        return;

    uint64_t const tail_selector_10 = (1 << 10) - 1;
    uint64_t const tail_selector_20 = (1 << 20) - 1;
    kmerID = code; // save trimmed kmer-code part
    uint64_t const k = (63 - __builtin_clzl(code)) >> 1;
    for (uint64_t i = 0; i < k - 4; ++i) // up-to four, i.e. 8 bits consecutive characters permitted
    {
        uint64_t tail_10 = tail_selector_10 & code;
        // XOR tail with A_4, C_4, G_4, T_
        if ((tail_10 == 0b00000000) || (tail_10 == 0b01010101) || (tail_10 == 0b10101010) || (tail_10 == 0b11111111))
        {
            // note that (x >> 64) won't be executed
            uint64_t offset = 64 - (std::max(PRIMER_MIN_LEN, k - i) - PRIMER_MIN_LEN);
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
                uint64_t offset = 64 - (std::max(PRIMER_MIN_LEN, k - i) - PRIMER_MIN_LEN);
                // delete all k bits ≥ len_selector
                prefix = (offset == 64) ? 0 : (prefix >> offset) << offset;
            }
        }
        if (!prefix)
            break;
        code >>= 2;
    }
    // add updated prefix
    kmerID += prefix;
}

// Check if not more than 3 out of the 5 last bases at the 3' end are CG.
// DNA sense/'+': 5' to 3' tests last 5 bps, antisense/'-': 3' to 5' tests first 5 bps
// Returns false if constraint is violated.
extern inline bool filter_CG_clamp(TKmerID const kmerID, char const sense, uint64_t const mask = 0)
{
    auto [prefix, code] = split(kmerID);
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

// primers should not end on TTT or ATT
extern inline bool filter_WWW_tail(TKmerID const kmerID, char const sense, uint64_t const mask = 0)
{
    if (sense != '+' && sense != '-')
        throw std::invalid_argument("Expected sense '-' or '+'!");
    std::unordered_set<uint64_t> const infixes = {0b000000, 0b000011, 0b001100, 0b001111, 0b110000, 0b110011, 0b111100, 0b111111};
    auto [prefix, code] = split(kmerID);
    uint64_t encoded_len = WORD_SIZE - __builtin_clzll(code) - 1;
    uint64_t tail;
    if (sense == '-') // reverse complement!
    {
        code >>= encoded_len - 6; // shift right to have prefix of length 10
        tail = code & 0b111111;
        return (infixes.find(tail) != infixes.end()) ? false : true;
    }
    if (mask)
    {
        uint64_t target_len = (PRIMER_MIN_LEN + __builtin_clzl(mask)) << 1;
        code >>= encoded_len - target_len; // delete prefix and trim to target length
    }
    return (infixes.find(code & 0b111111) != infixes.end()) ? false : true;
}

// Partial self-annealing test for leading PRIMER_MIN_LEN bps. Checks if any 4-mer
// occurs (reverse) complementary in first 16 bases. Case that suffix 4-mer is complementary
// to another suffix 4-mer is not tested yet.
// Deletes complete length mask if 4-more in reverse complement is found.
extern inline void filter_self_annealing(TKmerID & kmerID)
{
    // build 4-mer map of kmerID or kmerID2 if
    std::bitset<512> four;
    auto [prefix, code] = split(kmerID);
    uint8_t enc_l = (WORD_SIZE - __builtin_clzll(code) - 1);
    code >>= (enc_l - (PRIMER_MIN_LEN << 1));
    auto infix_mask = 0b11111111;
    while (code > (1 << 8))
    {
        four.set(infix_mask & code);
        code >>= 2;
    }
    code = kmerID & ~PREFIX_SELECTOR;
    // trim to last PRIMER_MIN_LEN bps
    // check annealing for same orientation
    uint8_t infix;
    while (code > (1UL << 8))
    {
        infix = ~(code & 0b11111111);
        if (four[infix])
        {
            uint8_t l = (WORD_SIZE - __builtin_clzll(code) - 1) >> 1;
            if (l <= PRIMER_MIN_LEN)
            {
                kmerID &= ~PREFIX_SELECTOR;
                return;
            }
            // del succeeding length bits because all longer kmers are affected
            kmerID = (kmerID & ~PREFIX_SELECTOR) + (prefix & ~((ONE_LSHIFT_63 >> (l - PRIMER_MIN_LEN - 1)) - 1));
            //std::cout << kmerID2str(kmerID) << std::endl;
        }
        code >>= 2;
    }
    // trim to last PRIMER_MIN_LEN bps
    // check reversed annealing
    code = kmerID & ~PREFIX_SELECTOR;
    int8_t offset = enc_l - 8;
    char const codes[4] = {'A', 'C', 'G', 'T'};
    uint64_t infixL;
    while (offset >= 0)
    {
        infixL = (0b11111111 << offset) & kmerID;
        infixL = 0b11111111 & ((~infixL) >> offset);
        infixL = ((infixL & 0b11) << 6) + ((infixL & 0b1100) << 2) + ((infixL & 0b110000) >> 2) + ((infixL & 0b11000000) >> 6);
        if (four[infixL])
        {
            uint8_t l = (WORD_SIZE - __builtin_clzll(0b11111111 << offset) - 1) >> 1;
            if (l <= PRIMER_MIN_LEN)
            {
                kmerID &= ~PREFIX_SELECTOR;
                return;
            }
            kmerID = (kmerID & ~PREFIX_SELECTOR) + (prefix & ~((ONE_LSHIFT_63 >> (l - PRIMER_MIN_LEN - 1)) - 1));
        }
        ----offset;
    }
}

// test forward and reverse cross annealing between two kmers. Two kmers cross-anneal
// if there are at least 4 nts complementary to each other.
extern inline void filter_cross_annealing(TKmerID & kmerID1, TKmerID & kmerID2)
{
    auto [prefix1, code1] = split(kmerID1);
    auto [prefix2, code2] = split(kmerID2);
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
    TKmerID code = kmerID & ~PREFIX_SELECTOR;
    assert(kmerID > 0 && code > 0);
    uint8_t AT = 0; // counter 'A'|'T'
    uint8_t CG = 0; // counter 'C'|'G'

    uint64_t enc_l;
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
        enc_l = (WORD_SIZE - 1 - __builtin_clzl(code)) >> 1;
        if ((code > (1ULL << 9)) && (((code & 0b11111111) == 0b11001100) || ((code & 0b11111111) == 0b00110011)))
        {
            if (enc_l <= PRIMER_MIN_LEN) // delete all length bits to discard this kmer
            {
                kmerID &= (1ULL << 52) - 1;
                return;
            }
            else // delete all length bits between current encoded length and subsequent
            {
                kmerID &= (ONE_LSHIFT_63 >> (enc_l - PRIMER_MIN_LEN + 1)) - 1;
            }
        }
        if (!(kmerID & PREFIX_SELECTOR))
            return;
        code >>= 2;
    }
    code = kmerID & ~PREFIX_SELECTOR;
    // start with largest kmer
    uint64_t tailing_zeros = ffsll((kmerID & PREFIX_SELECTOR) >> 54) - 1;
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
                {
                    kmerID ^= mask;
                }
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
    // Check for self-annealing
    filter_self_annealing(kmerID);
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

/* !\brief Check for self-dimerization, i.e. bonding energy by same sense bonding.
 * Returns true if ΔG ≥ -5kcal/mol
 */
extern inline bool filter_cross_dimerization(TKmerID kmerID1, TKmerID kmerID2)
{
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

// evaluate soft constraints
double score_kmer(TKmerID kmerID, uint64_t mask)
{
    double score = 0.0;
    double const eps = 0.1;
    // has 3 run
    auto [prefix, code] = split(kmerID);
    if (mask)
    {
        uint8_t enc_l = (WORD_SIZE - __builtin_clzll(code) - 1) >> 1;
        uint8_t target_l = __builtin_clzll(mask) + PRIMER_MIN_LEN;
        code >>= enc_l - target_l;
    }
    while (code >= (1 << 6))
    {
        uint64_t infix = code & 0b111111;
        if (infix == 0b000000 || infix == 0b111111)
            score -= eps;
        if (infix == 0b010101 || infix == 0b101010) // 3 runs of CG worse
            score -= 2*eps;
        code >>= 2;
    }
    return std::max(0.0, score);
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
