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

namespace priset
{

enum TIMEIT {
    MAP, // mappability computation
    TRANSFORM, //
    FILTER1, // chemical filter
    COMBINER, // kmer combiner
    FILTER2, // match chemical properties
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

//!\brief Enums for computational methods for primer melting temperature.
enum class TMeltMethod
{
    WALLACE,
    SALT_ADJUSTED
};

//using dna = typename seqan::Dna5;
typedef seqan::Dna5 dna;

// type declarations
using TSeqNo = uint64_t;
using TSeqPos = uint64_t;
// Kmer length type. A negative indicates reverse direction given associated position.
using TKmerLength = int64_t;
using TBWTLen = uint64_t;
using TFMIndexConfig = TGemMapFastFMIndexConfig<TBWTLen>;
typedef String<seqan::Dna, seqan::Alloc<>> TString;
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

/*
 * Datatype to store a kmer as an alphabet sequence, a melting temperature and
 * a unique integer ID.
 */
 // TODO: delete after transforming, only used in combiner
struct TKmer
{
    // Unique numerical identifier. Default 0 means unset.
    TKmerID ID{static_cast<TKmerID>(0)};

    // Melting temperature of k-mer sequence.
    float Tm{0};
};

// A bit vector for each reference, with 1 indicating a kmer starting position.
typedef sdsl::bit_vector TReference;
typedef std::vector<TReference> TReferences;

// Stores for each reference the encoded kmers in order of occurrence.
typedef std::vector<std::deque<TKmerID>> TKmerIDs;

// Transforms sequence ID 0..n-1 to contiguous range 0 .. k <= n-1.
// Background: some sequences may have no kmers and therefore no space should be reserved in bit vector set.
typedef std::unordered_map<TSeqNo, TSeqNo> TSeqNoMap;

// todo: TKmerLocation and TKmerPair inherit from same base struct.
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

 // Type for storing kmer combinations by their IDs and spatial occurences.
struct TKmerPair
{
    // The container type for storing a pair location, i.e. sequence ID, and start positions of fwd and rev primer.
    using TKmerPairLocation = typename std::tuple<TSeqNo, TSeqPos, TSeqPos>;
    using TKmerPairLocations = typename std::vector< TKmerPairLocation >;
    using const_iterator = TKmerPairLocations::const_iterator;
    using size_type = TKmerPairLocations::size_type;
private:
    // k-mer identifier for forward (5') primer sequence
    TKmerID kmer_ID1{0};
    // k-mer identifier for reverse (3') primer sequence
    TKmerID kmer_ID2{0};
    // absolute difference of their melting temperatures
    float Tm_delta;

public:
    // The container for storing pair locations.
    // TODO: make private and provide push_back and at functions
    TKmerPairLocations pair_locations;

    TKmerPair(TKmerID kmer_ID1_, TKmerID kmer_ID2_, float Tm_delta_) :
        kmer_ID1{kmer_ID1_}, kmer_ID2{kmer_ID2_}, Tm_delta{Tm_delta_} {}

    TKmerPair(TKmerID kmer_ID1_, TKmerID kmer_ID2_, float Tm_delta_, TKmerPairLocation & pair_location_) :
        kmer_ID1{kmer_ID1_}, kmer_ID2{kmer_ID2_}, Tm_delta{Tm_delta_}
    {
        pair_locations.push_back(pair_location_);
        // .resize(pair_locations_.size());
        //std::copy(pair_locations_.begin(), pair_locations.end(), pair_locations.begin());
    }

    TKmerID get_kmer_ID1() const
    {
        return kmer_ID1;
    }

    TKmerID get_kmer_ID2() const
    {
        return kmer_ID2;
    }

    float get_Tm_delta() const noexcept
    {
        return Tm_delta;
    }

    size_t container_size() const noexcept
    {
        return pair_locations.size();
    }

    TSeqNo accession_ID_at(size_t i) const noexcept
    {
        return std::get<0>(pair_locations[i]);
    }

    TSeqPos kmer_pos_at(size_type i, TKmerID kmerID = 1) const
    {
        if (kmerID == 1)
            return std::get<1>(pair_locations[i]);
        return std::get<2>(pair_locations[i]);
    }
};

// List type of pairs.
typedef std::vector<TKmerPair> TKmerPairs;

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

/*
* Struct to associate a taxonomic identifier with a set of accession numbers and sequences.
*/
template<typename taxid_type=unsigned int, typename accession_type=std::string,  typename sequence_type=std::string>
struct Reference
{
    //!\brief Taxonomic ID associated to the reference collection.
    const taxid_type taxid;

    //!\brief Default constructor.
    constexpr Reference() : taxid(static_cast<taxid_type>(1)) {}
    //!\brief Copy constructor.
    constexpr Reference(Reference const &) = default;
    //!\brief Copy construction via assignment.
    constexpr Reference & operator=(Reference const &) = default;
    //!\brief Move constructor.
    constexpr Reference (Reference &&) = default;
    //!\brief Move assignment.
    constexpr Reference & operator=(Reference &&) = default;
    //!\brief Use default deconstructor.
    ~Reference() = default;

    //!\brief Construct by taxid.
    explicit constexpr Reference(taxid_type const _taxid) noexcept : taxid(_taxid) {}

    void insert(accession_type const accession, sequence_type const sequence)
    {
        assert(accession.size() > 0 && sequence.size() > 0);
        _accessions.push_back(accession);
        _sequences.push_back(sequence);
    }

    //!\brief Erase a sequence given its accession.
    bool erase(accession_type accession)
    {
        auto it = std::find(_accessions.begin(), _accessions.end(), accession);
        if (it != _accessions.end())
        {
            size_t offset = it - _accessions.begin();
            _accessions.erase(it);
            _sequences.erase(_sequences.begin() + offset);
        }
        else
            return false;
        return true;
    }

private:
    //!\brief Accessions associated to taxid.
    std::vector<accession_type> _accessions;
    //!\brief List of sequences associated to the accessions.
    std::vector<sequence_type> _sequences;
};

} // priset
