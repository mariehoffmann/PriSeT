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
#include <vector>

#include <seqan/basic.h>

//#include "seqan::Dna5.hpp"
//#include "primer_config.hpp"
#include "types.hpp"

// satisfies the primer_config_concept.
namespace priset::chemistry
{

//!\brief Wallace rule to compute the melting temperature of a primer sequence.
float primer_melt_wallace(seqan::String<priset::dna> const & sequence)
{
    assert(length(sequence) < (1 << 8));
    uint8_t cnt_AT = 0;
    uint8_t CG_cnt = 0;
    for (auto c : sequence)
    {
        switch(char(c))
        {
            case 'A':
            case 'T': ++cnt_AT; break;
            case 'C':
            case 'G': ++CG_cnt; break;
            default: std::cout << "ERROR: unsupported sequence character '" << c << "'\n";
        }
    }
    return 2*cnt_AT + 4*CG_cnt;
}

//!\brief Salt-adjusted method to compute melting temperature of primer sequence.
// input primer:string sequence, Na:float molar Natrium ion concentration
inline float primer_melt_salt(seqan::String<priset::dna> const & sequence, float const Na)
{
    assert(length(sequence) < (1 << 8));
    uint8_t cnt_AT = 0;
    uint8_t CG_cnt = 0;
    for (auto c : sequence)
    {
        switch(char(c))
        {
            case 'A':
            case 'T': ++cnt_AT; break;
            case 'C':
            case 'G': ++CG_cnt; break;
            default: std::cout << "ERROR: unsupported sequence character '" << c << "'\n";
        }
    }
    return 100.5 + 41.0 * CG_cnt / seqan::length(sequence) - 820.0 / seqan::length(sequence) + 16.6 * std::log10(Na);
}

//!\brief Compute melting temperature of primer sequence with method set in primer configuration.
template<typename primer_config_type>
float get_Tm(primer_config_type const & primer_cfg, seqan::String<priset::dna> const & sequence) noexcept
{
    switch(primer_cfg.get_primer_melt_method())
    {
        case WALLACE: return primer_melt_wallace(sequence);
        default: return primer_melt_salt(sequence, primer_cfg.get_Na());
    }
}

//!\brief Check if melting temperature is in range set by the primer configurator.
template<typename primer_config_type>
inline bool filter_Tm(primer_config_type const & primer_cfg, seqan::String<priset::dna> const & sequence)
{
    float Tm = get_Tm(primer_cfg, sequence);
    if (Tm >= primer_cfg.get_min_Tm() && Tm <= primer_cfg.get_max_Tm())
        return true;
    return false;
}

//!\brief Check if CG content is in range set by the primer configurator.
// Returns false if constraint is violated.
template<typename primer_cfg_type>
inline bool filter_CG(primer_cfg_type const & primer_cfg, seqan::String<priset::dna> const & sequence)
{
    assert(length(sequence) < (1 << 8));
    uint8_t CG_cnt = 0;
    for (auto c : sequence)
    {
        switch(char(c))
        {
            case 'A': break;
            case 'C':
            case 'G': ++CG_cnt; break;
            case 'T': break;
            default: std::cout << "ERROR: unsupported sequence character '" << c << "'\n";
        }
    }
    float CG = float(CG_cnt) / float(seqan::length(sequence));
    return (CG >= primer_cfg.get_min_CG_content() && CG <= primer_cfg.get_max_CG_content());
}

//!\brief Check if not more than 3 out of the 5 last bases at the 3' end are CG.
//  DNA sense/'+': 5' to 3', antisense/'-': 3' to 5'
// Returns false if constraint is violated.
template<typename primer_config_type>
inline bool filter_CG_clamp(primer_config_type const & primer_cfg, seqan::String<priset::dna> const & sequence, char const sense)
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
 * Returns true if generation of secondary structures is improbable.
 * Tested are hairpins, self- and cross dimerization.
 */
template<typename primer_cfg_type>
inline bool filter_secondary_structures(primer_cfg_type const & primer_cfg)
{
    return true;
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

/*
# one-letter encoding for set of aligned sequences, no gaps
def compress_helper(aligned_sequences, pos, length, bin_codes):
    seq_x = ''
    for i in range(pos, pos+length):
        code = reduce(lambda x, y: x|y, [bin_codes[aseq.seq[i]] for aseq in aligned_sequences])
        seq_x += one_letter_encode[code]
    return seq_x

def compress(aligned_sequences, pos, length):
    return compress_helper(aligned_sequences, pos, length, {'A': 1, 'C': 2, 'G': 4, 'T': 8})

def complement_compress(aligned_sequences, pos, length):
    return compress_helper(aligned_sequences, pos, length, {'A': 8, 'C': 4, 'G': 2, 'T': 1})[::-1]
*/

} // namespace chemistry
