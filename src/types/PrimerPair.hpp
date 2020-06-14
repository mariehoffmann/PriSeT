#pragma once

#include "../dna.hpp"
#include "CombinePattern.hpp"
#include "simple_types.hpp"

namespace priset
{

/*
 * A pair of encoded kmers (kmer IDs) given by their indices in the kmer ID
 * container associated to a reference bit vector. Since kmer IDs may encode up
 * to |[primer_max_length : primer_min_length]| kmers, the combine pattern stores
 * which kmers are combined. See encoding scheme here: CombinePattern.
 */
struct PrimerPair
{
    PrimerPair() = default;
    // Constructor with parameters. Note that a valid primer pair cannot be built
    // from two k-mers starting at the same sequence position, i.e. they have the
    // same rank in the common reference!
    PrimerPair(uint64_t _seqNo_cx, uint64_t _r_fwd, uint64_t _r_rev, CombinePattern _cp) :
        seqNo_cx(_seqNo_cx), r_fwd(_r_fwd), r_rev(_r_rev), cp(_cp)
    {
            assert(r_fwd < r_rev);
    }

    ~PrimerPair() = default;

    // Get compressed sequence identifier where this pair originates from.
    TSeqNo get_seqNo() const
    {
        // ranks are uninitialized if both zeros
        assert(r_fwd || r_rev);
        return seqNo_cx;
    }

    // Get rank for forward TKMerID.
    uint64_t get_rank_fwd() const
    {
        // ranks are uninitialized if both zeros
        assert(r_fwd || r_rev);
        return r_fwd;
    }

    // Get rank for forward TKMerID.
    uint64_t get_rank_rev() const
    {
        // ranks are uninitialized if both zeros
        assert(r_fwd || r_rev);
        return r_rev;
    }

    CombinePattern get_combine_pattern() const noexcept
    {
        return cp;
    }

private:
    // Sequence identifier (equivalent to sequence position in text corpus).
    TSeqNo seqNo_cx{0};

    // Rank of forward kmer.
    uint64_t r_fwd{0};

    // Rank of reverse kmer.
    uint64_t r_rev{0};

    // Length combinations of kmers as bit mask.
    CombinePattern cp{};
};

}  // namespace priset
