#pragma once

#include <vector>

#include "GenMapTypes.hpp"

namespace priset
{

// vector type of k-mers and their locations
struct KmerLocation
{
private:
    // The unique kmer identifier encoding dna4 sequence up to 30 bp.
    TKmerID kmer_ID{0};
    // TODO: delete K
    TKmerLength K{0};

public:
    using TLocationVec = typename std::vector<TLocation>;

    TLocationVec locations{}; // TODO: make private and provide public const_iterator for it

    KmerLocation(TKmerID kmer_ID_, TKmerLength K_, TLocationVec & locations_) :
    kmer_ID{kmer_ID_}, K{K_}
    {
        locations.resize(locations_.size());
        std::copy(locations_.begin(), locations_.end(), locations.begin());
    }

    using size_type = TLocationVec::size_type;
    using const_iterator = TLocationVec::const_iterator;

    void set_kmer_ID(TKmerID kmer_ID_) noexcept
    {
        kmer_ID = kmer_ID_;
    }

    TKmerID get_kmer_ID() const noexcept
    {
        return kmer_ID;
    }

    TKmerID get_kmer_ID1() const noexcept
    {
        return get_kmer_ID();
    }

    TKmerID get_kmer_ID2() const noexcept
    {
        return 0;
    }

    TKmerLength get_K() const noexcept
    {
        return K;
    }

    size_type container_size() const noexcept
    {
         return locations.size();
    }

    TSeqNo accession_ID_at(size_type i) const
    {
         return seqan::getValueI1<TSeqNo, TSeqPos>(locations[i]);
    }

    TSeqPos kmer_pos_at(size_type i) const
    {
        return seqan::getValueI2<TSeqNo, TSeqPos>(locations[i]);
    }

};

}  // namespace priset
