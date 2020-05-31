#pragma once

#include <vector>

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
    Pair(uint64_t reference_, uint64_t r_fwd_, uint64_t r_rev_, CombinePattern cp_) :
        reference(reference_), r_fwd(r_fwd_), r_rev(r_rev_), cp(cp_) {}
    //std::tuple<uint64_t, uint64_t, CombinePattern<TKmerID, TKmerLength>> Pair;
    ~Pair() = default;

    // Reference identifier (equivalent to position in corpus).
    uint64_t reference;
    // Rank of forward kmer.
    uint64_t r_fwd;
    // Rank of reverse kmer.
    uint64_t r_rev;
    // Length combinations of kmers as bit mask.
    CombinePattern cp;
};

template<typename Pair>
using PairList = std::vector<Pair>;

}

struct PairUnpacked
{
    PairUnpacked() = default;
    PairUnpacked(size_t const _reference_count) : reference_count(_reference_count)
    {
        references.resize(reference_count, 0);
    };

    // Set flag in reference vector if this pair has a match in a reference
    // with ID = seqNo_cx.
    void set_reference_match(size_t const seqNo_cx)
    {
        if (references.size() <= seqNo_cx)
            references.resize(seqNo_cx + 1);
        references[seqNo_cx] = 1;
    }

    // Get flag for match of reference sequence with ID = seqNo_cx.
    bool get_reference_match(size_t seqNo_cx) const noexcept
    {
        return ((references.size() <= seqNo_cx) ? 0 : references[seqNo_cx];
    }

    // Return
    size_t get_reference_count()

private:
    TKmerID fwd;
    uint64_t mask_fwd;
    TKmerID rev;
    uint64_t mask_rev;

    size_t reference_count;
    std::vector<bool> references;
};

}
