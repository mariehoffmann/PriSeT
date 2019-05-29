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

// map of k-mer locations (determined by genmap)
typedef std::map<seqan::Pair<priset::TSeqNo, priset::TSeqPos>,
         std::pair<std::vector<seqan::Pair<priset::TSeqNo, priset::TSeqPos> >,
                   std::vector<seqan::Pair<priset::TSeqNo, priset::TSeqPos> > > > TLocations;

//
using TDirectoryInformation = typename seqan::StringSet<seqan::CharString, seqan::Owner<seqan::ConcatDirect<> > > ;
// container for fasta header lines
using TSequenceNames = typename seqan::StringSet<seqan::CharString, seqan::Owner<seqan::ConcatDirect<> > >;
// container for fasta sequence lengths
using TSequenceLengths = typename seqan::StringSet<uint32_t>;

typedef uint64_t TKmerID;
/*
 * Datatype to store a kmer as an alphabet sequence, a melting temperature and
 * a unique integer ID.
 */
struct TKmer
{
    // DNA sequence it represents
    TSeq seq{};

    // Unique numerical identifier.
    TKmerID ID{static_cast<TKmerID>(1<<24)};

    // Melting temperature of k-mer sequence.
    float Tm{0};
};

// store kmer combinations by their IDs
typedef std::vector<std::pair< TKmerID, TKmerID > > TPairs;

// vector type of k-mers and their locations
typedef std::vector<std::pair<TKmer, std::vector<seqan::Pair<priset::TSeqNo, priset::TSeqPos> > > > TKmerLocations;

/*
 * Datatype to store matches of two k-mers within a list of accessions.
 */
template<typename kmer_type>
struct match
{
private:

    // k-mer identifier for forward (5') primer sequence
    kmer_type kmer_fwd;
    // k-mer identifier for reverse (3') primer sequence
    kmer_type kmer_rev;
    // taxid of all accessions in the accession list
    typename kmer_type::primer_cfg_type::taxid_type taxid;
    // absolute difference of their melting temperatures
    typename kmer_type::primer_cfg_type::float_type Tm_delta;
    // accessions from library where both k-mers match and which are directly assigned
    // to the taxid (e.g. taxid is not the lca)
    // TODO: decide to store here also the location asscociated with an accession
    std::vector<typename kmer_type::primer_cfg::accession_type> accession_list;
public:
    //using = ;
    constexpr match() = default;
    match(kmer_type kmer_fwd, kmer_type kmer_rev)
    {
        Tm_delta = std::abs(kmer_fwd.get_Tm() - kmer_fwd.getTm());
    }
    ~match() = default;
};

// TODO: use key value tuple and use ordered set or map
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
