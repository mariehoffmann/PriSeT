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
