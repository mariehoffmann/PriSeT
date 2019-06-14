// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================
//          Author: Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
//          Manual: https://github.com/mariehoffmann/PriSeT

// Interval type for storing upper and lower bounds.

#pragma once

#include <bitset>
#include <stdlib.h>

#include "../submodules/genmap/src/common.hpp"

//#include "primer_cfg_type.hpp"

namespace priset
{

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

//
using TDirectoryInformation = typename seqan::StringSet<seqan::CharString, seqan::Owner<seqan::ConcatDirect<> > > ;

// container for fasta header lines
using TSequenceNames = typename seqan::StringSet<seqan::CharString, seqan::Owner<seqan::ConcatDirect<> > >;

// container for fasta sequence lengths
using TSequenceLengths = typename seqan::StringSet<uint32_t>;

// The type for kmer identifiers, 1-based.
typedef uint64_t TKmerID;

// The type for taxonomic identifiers.
typedef uint32_t TTaxid;

// The type for numerical accession identifiers, 1-based.
typedef uint64_t TAccID;

/*
 * Datatype to store a kmer as an alphabet sequence, a melting temperature and
 * a unique integer ID.
 */
struct TKmer
{
    // Unique numerical identifier. Default 0 means unset.
    TKmerID ID{static_cast<TKmerID>(0)};

    // alphabet sequence of k-mer
    TSeq seq{};

    // Melting temperature of k-mer sequence.
    float Tm{0};
};

// The map to resolve kmer IDs and their structs.
// TODO: ID redundant, see application to possibly remove from struct
typedef std::map<TKmerID, TKmer> TKmerMap;

// vector type of k-mers and their locations
typedef std::pair<TKmerID, std::vector<TLocation > > TKmerLocation;
typedef std::vector<TKmerLocation > TKmerLocations;

 // Type for storing kmer combinations by their IDs and spatial occurences.
struct TPair
{
    using TPairLocations = typename std::vector<std::tuple<TSeqNo, TSeqPos, TSeqPos> >;
    // k-mer identifier for forward (5') primer sequence
    TKmerID kmer_fwd;
    // k-mer identifier for reverse (3') primer sequence
    TKmerID kmer_rev;
    // absolute difference of their melting temperatures
    float Tm_delta;
    // The set of locations given by sequence id and position indices of fwd and rev sequence IDs.
    TPairLocations pair_locations;
};

// List type of pairs.
typedef std::vector<TPair> TKmerPairs;

// Result table output for app
struct TResult
{
    // taxonomic node identifier
    TTaxid taxid;
    // forward kmer identifier
    TKmerID fwd;
    // reverse kmer identifier
    TKmerID rev;
    // taxonomic node counter with suitable fwd/rev primer matches in their descendents
    uint16_t match_ctr;
    // total number of taxonomic nodes in subtree having non-zero accessions assigned
    uint16_t covered_taxids;
    // list of references (empty for taxids having no direct assignments)
    std::vector<std::string> accs;
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
