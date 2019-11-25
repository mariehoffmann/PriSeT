// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================
//          Author: Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
//          Manual: https://github.com/mariehoffmann/PriSeT

// Interval type for storing upper and lower bounds.

#pragma once

#include <bitset>
#include <stdlib.h>

#include <sdsl/bit_vectors.hpp>

#include "../submodules/genmap/src/common.hpp"

//#include "primer_cfg_type.hpp"

#define ONE_LSHIFT_63 9223372036854775808ULL

#define NULL_KMERID 0ULL

namespace priset
{

enum TIMEIT {
    MAP, // mappability computation
    FILTER1_TRANSFORM, //
    COMBINE_FILTER2, // kmer combiner
    PAIR_FREQ, // Pair frequency cutoff
    SIZE
};

enum KMER_COUNTS
{
    MAP_CNT, // single kmers
    FILTER1_CNT, // kmer pairs
    COMBINER_CNT,
    FILTER2_CNT
};

typedef std::array<uint64_t, 4> TKmerCounts;

//using dna = typename seqan::Dna5;
typedef seqan::Dna5 dna;

// type declarations
using TSeqNo = uint64_t;
using TSeqPos = uint64_t;
// Kmer length type. A negative indicates reverse direction given associated position.
using TKmerLength = int64_t;
using TBWTLen = uint64_t;
using TFMIndexConfig = TGenMapFastFMIndexConfig<TBWTLen>;
typedef seqan::String<seqan::Dna, seqan::Alloc<>> TString;
typedef seqan::StringSet<TString, seqan::Owner<seqan::ConcatDirect<SizeSpec_<TSeqNo, TSeqPos> > > > TStringSet;
// set index type, TBiIndexConfig defined src/common.hpp
using TIndex = seqan::Index<TStringSet, TBiIndexConfig<TFMIndexConfig> >;

typedef seqan::String<priset::dna> TSeq;

// The location type defined by sequence ID and position.
typedef seqan::Pair<priset::TSeqNo, priset::TSeqPos> TLocation;

// The map of k-mer locations.
typedef std::map<TLocation,
         std::pair<std::vector<TLocation >,
                   std::vector<TLocation > > > TLocations;

// A k-mer location augmented by the information about K.
typedef std::tuple<priset::TSeqNo, priset::TSeqPos, priset::TKmerLength> TKLocation;

// The map of k-mer locations augmented by K information to preserve key uniqueness.
// Contains mappings for all values of K in range (see primer_cfg_type.hpp).
typedef std::map<TKLocation,
        std::pair<std::vector<TLocation >,
                  std::vector<TLocation > > > TKLocations;

//
using TDirectoryInformation = typename seqan::StringSet<seqan::CharString, seqan::Owner<seqan::ConcatDirect<> > > ;

// container for fasta header lines
using TSequenceNames = typename seqan::StringSet<seqan::CharString, seqan::Owner<seqan::ConcatDirect<> > >;

// container for fasta sequence lengths
using TSequenceLengths = typename seqan::StringSet<uint32_t>;

// The type for kmer identifiers encoding kmer sequences up to a length of 30 bp.
typedef uint64_t TKmerID;

// The type for taxonomic identifiers.
typedef uint32_t TTaxid;

// The type for numerical accession identifiers, 1-based.
typedef uint64_t TAccID;

// The type of an accession
typedef std::string TAcc;

// A bit vector for each reference, with 1 indicating a kmer starting position.
typedef sdsl::bit_vector TReference;
typedef std::vector<TReference> TReferences;

// Stores for each reference the encoded kmers in order of occurrence.
typedef std::vector<std::deque<TKmerID>> TKmerIDs;

// Transforms sequence ID 0..n-1 to contiguous range 0 .. k <= n-1.
// Background: some sequences may have no kmers and therefore no space should be reserved in bit vector set.
typedef std::unordered_map<TSeqNo, TSeqNo> TSeqNoMap;

// vector type of k-mers and their locations
struct TKmerLocation
{
private:
    // The unique kmer identifier encoding dna4 sequence up to 30 bp.
    TKmerID kmer_ID{0};
    // TODO: delete K
    TKmerLength K{0};

public:
    using TLocationVec = typename std::vector<TLocation>;

    TLocationVec locations{}; // TODO: make private and provide public const_iterator for it

    TKmerLocation(TKmerID kmer_ID_, TKmerLength K_, TLocationVec & locations_) :
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

struct TPrimerPair
{
    std::string name_fwd;
    uint64_t code_fwd;
    std::string name_rev;
    uint64_t code_rev;
};

// Result table output for app
struct TResult
{
    // taxonomic node identifier
    TTaxid taxid;
    // forward kmer identifier
    TKmerID kmer_ID1;
    // reverse kmer identifier
    TKmerID kmer_ID2;
    // taxonomic node counter with suitable fwd/rev primer matches in their descendents
    uint16_t match_ctr;
    // total number of taxonomic nodes in subtree having non-zero accessions assigned
    uint16_t covered_taxids;
    // list of references (empty for taxids having no direct assignments)
    std::vector<std::string> accIDs{};
    // Return comma-separated result string
    // Q: can R read csv tables with varying columns? if not use array notation
    std::string to_string()
    {
        std::stringstream ss;
        ss << taxid << "," << kmer_ID1 << "," << kmer_ID2 << "," << match_ctr << "," << covered_taxids;
        if (accIDs.size())
            for (auto accID : accIDs)
                ss << "," << accID;
        ss << "\n";
        return ss.str();
    }
};

// Upstream result collection as map. Since tuples are not hashable, it is converted into a string before hashing.
struct TUpstreamKey
{
    using THash = std::string;
    TUpstreamKey(TTaxid taxid_, TKmerID fwd_, TKmerID rev_) : taxid(taxid_), fwd(fwd_), rev(rev_)
    {
    }

    // String represenation, usable as dictionary key or direct csv output.
    THash to_string()
    {
        std::stringstream s;
        s << taxid << "," << fwd << "," << rev;
        return s.str();
    }

    TTaxid taxid;
    TKmerID fwd;
    TKmerID rev;
};

} // priset
