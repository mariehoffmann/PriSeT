// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================

/*!\file
 * \brief Global Functions for Primer Constraint Checking.
 * \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <cmath>

#include <priset/core/PrimerConfig>
#include <priset/types/dna.hpp>

// satisfies the primer_config_concept.
namespace priset::chemistry
{
//!\brief Enums for computational methods for primer melting temperature.
enum method {wallace, salt_adjusted};

//!\brief Wallace rule to compute the melting temperature of a primer sequence.
template<typename sequence_type, typename float_type>
// todo: require base type of string compatible with char via requires concept
//sequence_type::value_type == dna
float_type primer_melt_wallace(sequence_type primer, sequence_type::iterator it1, sequence_type::iterator it2)
{
    assert((std::is_same<sequence_type::value_type, dna>));
    size_t cnt_AT = std::count_if(it1, it2, [](dna b) {return b == dna::A || c == dna::T;});
    size_t cnt_CG = std::count_if(it1, it2, [](dna b) {return c == dna::C || c == dna::G;});
    return 2*cnt_AT + 4*cnt_CG;
}

//!\brief Salt-adjusted method to compute melting temperature of primer sequence.
// input primer:string sequence, Na:float molar Natrium ion concentration
template<typename sequence_type, typename float_type>
float_type primer_melt_salt(sequence_type primer, float_type Na, sequence_type::iterator it1, sequence_type::iterator it2)
{
    float_type cnt_CG = std::count_if(it1, it2, \
        [](char c) {return c == 'C' || c == 'G';});
    sequence_type::size_type primer_len = static_cast<sequence_type::size_type>(it2-it1);
    return 100.5 + 41.0*cnt_CG/primer_len - 820.0/primer_len + 16.6*std::log10(Na);
}

//!\brief  variation measure for a block of aligned sequences in terms of number of columns
// having different nucleotides. A low score indicates high conservation.
// input sequences:[[]], pos:int, offset:int
/*def variation_score(aligned_sequences, pos, offset):
    transposed = [[aseq.seq[i] for aseq in aligned_sequences] for i in range(pos, pos+offset)]
    score = sum([1 if len(set(col)) > 1 else 0 for col in transposed])
    return score
*/

//!\brief Compress a window of aligned sequences to 1-letter encode.
template<typename sequence_type>
// todo: use aligned sequence type
// block needs to be gap-free, N is ignored due to its ambiguity, all sequences are padded to the same length
std::vector<dna> block_compress(std::vector<sequence_type> const aligned_sequences,
    size_t const pos, size_t const offset)
{
    // assert that block end doesn't exceed the aligned sequences lengths
    assert(aligned_sequences[0].size() - pos >= offset);
    //matchstr = ['N' for _ in range(min(offset, len(aligned_sequences[0].seq)-pos))]
    std::vector<dna> as_cx(offset);
    //codes = {1: '|', 2: '2', 3: '3', 4: '4'}
    unsigned short int mask_A = 1, mask_C = 2, mask_G = 4, mask_T = 8;
    for (unsigned int i = offset; i < offset + size; ++i)
    {
        for (unsigned int j = 0; j < aligned_sequences.size(); ++j)
        {
            unsigned short int column_mask = 0;
            switch(aligned_sequences[i][j])
            {
                case 'A': column_mask |= mask_A; break;
                case 'C': column_mask |= mask_C; break;
                case 'G': column_mask |= mask_G; break;
                case 'T': column_mask |= mask_T; break;
                default: std::cout << "Warning: block_compress scans block with unknown symbol '" <<
                    aligned_sequences[i][j] << "'" << std::endl;
            }
        }
        // build dna string conform to key set of priset::str2dna map
        std::string dna_str = "";
        if (column_mask & mask_A == mask_A)
            dna_str.append('A');
        if (column_mask & mask_C == mask_C)
            dna_str.append('C');
        if (column_mask & mask_G == mask_G)
            dna_str.append('G');
        if (column_mask & mask_T == mask_T)
            dna_str.append('T');

        as_cx[i-offset] = str2dna[dna_str];
    }
    return as_cx
}

//!\brief Computer melting temperature of primer sequence.
template<typename sequence_type, typename float_type>
float_type get_Tm(sequence_type const primer, PrimerConfig& const primer_cfg) noexcept
{
    switch(method primer_cfg.get_primer_melt_method())
    {
        case wallace: return primer_melt_wallace(primer);
        default: return primer_melt_salt(primer, primer_cfg.get_Na());
    }
}

//!\brief Check for all target sequences the melting temperature range.
template<typename sequence_type, typename float_type>
bool filter_Tm(PrimerConfig& const primer_cfg, sequence_type& const sequence)
{
    return filter_Tm(primer_cfg, sequence, 0, sequence.size());
}

//!\brief Check melting temperature for subsequence.
template<typename sequence_type, typename float_type>
bool filter_Tm(PrimerConfig& const primer_cfg, sequence_type& const sequence,
    sequence_type::size_type const offset, sequence_type::size_type const size)
{
    float_type Tm = get_Tm(primer, primer_cfg, offset, size);
    if (Tm >= primer_cfg.get_min_Tm() && Tm <= primer_cfg.get_max_Tm())
        return true;
    return false;
}

//!\brief Compute relative GC content.
// TODO: require sequence_concept (includes sequence has iterator)
template<typename sequence_type, typename float_type>
float_type get_GC_content(sequence_type sequence, PrimerConfig & cfg,
    sequence_type::iterator_type it1, sequence_type::iterator_type it2) noexcept
{
    size_t cnt_GC = std::count_if(it1, it2, [](dna b) {return b == dna::C || b == dna::G;});
    return static_cast<float_type>(cnt_GC)/static_cast<float_type>(it2-it1);
}

//!\brief Check if GC content is in the recommended range.
template<typename sequence_type, typename float_type>
bool filter_GC_content(sequence_type sequence, PrimerConfig & cfg,
    sequence_type::size_type offset, sequence_type::size_type size) noexcept
{
    float_type CG_content = get_CG_content(sequence, cfg, offset, size);
    if (CG_content < cfg.get_min_CG_content() || CG_content > cfg.get_max_CG_content())
        return false;
    return true;
}

/*
# check for GC at 3' end, DNA sense/'+': 5' to 3', antisense/'-': 3' to 5', should be <= 3 in last 5 bps
def filter_GC_clamp(sequence, sense='+'):
    if sense == '+':
        gc = len([1 for nt in sequence[-5:] if nt in ['C', 'G']])
    else:
        gc = len([1 for nt in sequence[:5] if nt in ['C', 'G']])
    if gc > 3:
        return False, gc
    return True, gc

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

} // namespace priset::chemistry
