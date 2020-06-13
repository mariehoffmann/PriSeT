#pragma once

#include "../dna.hpp"
#include "simple_types.hpp"

namespace priset
{

/*
 * A pair of encoded kmers (kmer IDs) given by their indices in the kmer ID
 * container associated to a reference bit vector. Since kmer IDs may encode up
 * to |[primer_max_length : primer_min_length]| kmers, the combine pattern stores
 * which kmers are combined. See encoding scheme here: CombinePattern.
 */
template<typename CombinePattern>
struct PrimerPair
{
    PrimerPair() = default;
    PrimerPair(uint64_t _seqNo_cx, uint64_t _r_fwd, uint64_t _r_rev, CombinePattern _cp) :
        seqNo_cx(_seqNo_cx), r_fwd(_r_fwd), r_rev(_r_rev), cp(_cp) {}

    ~PrimerPair() = default;

    // Get compressed sequence identifier where this pair originates from.
    TSeqNo get_seqNo() const
    {
        return seqNo_cx;
    }

    // Get rank for forward TKMerID.
    uint64_t get_rank_fwd() const
    {
        return r_fwd;
    }

    // Get rank for forward TKMerID.
    uint64_t get_rank_rev() const
    {
        return r_rev;
    }

    CombinePattern get_combine_pattern() const
    {
        return cp;
    }

private:
    // Sequence identifier (equivalent to sequence position in text corpus).
    TSeqNo seqNo_cx;

    // Rank of forward kmer.
    uint64_t r_fwd;

    // Rank of reverse kmer.
    uint64_t r_rev;

    // Length combinations of kmers as bit mask.
    CombinePattern cp;
};

}  // namespace priset
