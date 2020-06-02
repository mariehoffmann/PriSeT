#pragma once

#include <vector>

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
struct Pair
{
    Pair() = default;
    Pair(uint64_t _reference_ID, uint64_t _r_fwd, uint64_t _r_rev, CombinePattern _cp) :
        reference_ID(_reference_ID), r_fwd(_r_fwd), r_rev(_r_rev), cp(_cp) {}
    //std::tuple<uint64_t, uint64_t, CombinePattern<TKmerID, TKmerLength>> Pair;
    ~Pair() = default;

    // Reference identifier (equivalent to position in corpus).
    uint64_t reference_ID;

    // Rank of forward kmer.
    uint64_t r_fwd;

    // Rank of reverse kmer.
    uint64_t r_rev;

    // Length combinations of kmers as bit mask.
    CombinePattern cp;
};

template<typename Pair>
using PairList = std::vector<Pair>;

struct PairUnpacked
{
    PairUnpacked() = default;
    PairUnpacked(size_t const _reference_count) : reference_count(_reference_count)
    {
        reference_IDs.resize(reference_count, 0);
    };

    // Set flag in reference vector if this pair has a match in a reference
    // with ID = seqNo_cx.
    void set_reference_match(uint64_t const seqNo_cx)
    {
        if (reference_IDs.size() <= seqNo_cx)
            reference_IDs.resize(seqNo_cx + 1);
        reference_IDs[seqNo_cx] = 1;
    }

    // Get flag for match of reference sequence with ID = seqNo_cx.
    bool get_reference_match(size_t seqNo_cx) const noexcept
    {
        return (reference_IDs.size() <= seqNo_cx) ? 0 : reference_IDs[seqNo_cx];
    }

    // Return number of distinct reference_IDs
    size_t get_reference_count() noexcept
    {
        if (reference_count)
            return reference_count;
        reference_count = std::accumulate(reference_IDs.begin(), reference_IDs.end(), 1);
        return reference_count;
    }

private:
    TKmerID fwd;
    uint64_t mask_fwd;
    TKmerID rev;
    uint64_t mask_rev;

    size_t reference_count{0};
    std::vector<bool> reference_IDs;
};

}
