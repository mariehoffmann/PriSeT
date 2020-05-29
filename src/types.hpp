// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================
//          Author: Marie Hoffmann <ozymandiaz147 AT gmail.com>
//          Manual: https://github.com/mariehoffmann/PriSeT

// Interval type for storing upper and lower bounds.

#pragma once

#include <bitset>
#include <stdlib.h>
#include <vector>

#include "utilities.hpp"

#include "../submodules/sdsl-lite/include/sdsl/bit_vectors.hpp"
#include "../submodules/genmap/src/common.hpp"

#define ONE_LSHIFT_63 9223372036854775808ULL

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

static TKmerID NULL_TKMERID = 0;

// The type for taxonomic identifiers.
typedef uint32_t uint64_t;

// The type for numerical accession identifiers, 1-based.
typedef uint64_t TAccID;

// The type of an accession
typedef std::string TAcc;

// A bit vector for each reference, with 1 indicating a kmer starting position.
typedef sdsl::bit_vector TReference;
typedef std::vector<TReference> TReferences;

// Stores for each reference the encoded kmers in order of occurrence.
typedef std::vector<std::deque<TKmerID>> TKmerIDs;

// Translates sequences identifiers (seqNo) in use to a contiguous range (seqNo_cx).
// Dictionary is bidirectional: seqNo -> seqNo_cx and inverse add a leading one
// to the compressed key: (1 << 63 | seqNo_cx) -> seqNo.
// Background: some sequences produce no k-mers and therefore no space should be
// reserved in its bit transformation.
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


public:
    TPrimerPair(uint64_t _code_fwd, uint64_t _code_rev) : code_fwd(_code_fwd), code_rev(_code_rev){}

    uint64_t get_code_fwd() const
    {
        return code_fwd;
    }

    uint64_t get_code_rev() const
    {
        return code_rev;
    }

    std::string get_seq_fwd() const
    {
        return dna_decoder(code_fwd);
    }

    std::string get_seq_rev() const
    {
        return dna_decoder(code_rev);
    }

    std::pair<float, float> get_Tm()
    {
        return std::pair<float, float>{};
    }

};

// Result type storing a single primer pair with additional information.
struct TResult
{

private:

    // List of distinct taxonomic node identifiers in which this pair occurs.
    std::vector<uint64_t> taxa;

    // List of distinct references in which this pair occurs.
    std::vector<uint64_t> references;

    std::string primer_fwd;
    std::string primer_rev;

    int Tm_fwd{0};
    int Tm_rev{0};

    float CG_fwd{0};
    float CG_rev{0};

public:
    TResult(TKmerID kmerID_fwd, uint64_t mask_fwd, TKmerID _kmerID_rev, uint64_t _mask_rev, std::vector<uint64_t> _taxa, std::vector<uint64_t> _references) :
    taxa(_taxa), references(_references)
    {
        primer_fwd = dna_decode(kmerID_fwd, mask_fwd);
        primer_rev = dna_decode(kmerID_rev, mask_rev);
        Tm_fwd = get_Tm(kmerID_fwd, mask_fwd);
        Tm_rev = get_Tm(kmerID_rev, mask_rev);
        CG_fwd = CG(kmerID_fwd, mask_fwd);
        CG_rev = CG(kmerID_rev, mask_rev);
    }

    std::pair<string, string> get_primer_names() const
    {
        return std::pair<string, string>{hash<string>{}(primer_fwd), hash<string>{}(primer_rev)};
    }

    std::pair<string, string> get_primer_sequences() const
    {
        return std::pair<string, string>{primer_fwd, primer_rev};
    }

    std::vector<uint64_t> get_taxa() const
    {
        return taxa;
    }

    std::vector<uint64_t> get_references() const
    {
        return references;
    }

    // Return comma-separated result string
    std::string to_string()
    {
        std::stringstream ss;
        ss << primer_fwd << "," << Tm_fwd << "," << CG_fwd << "," << primer_rev << "," << Tm_rev << "," << CG_rev << ",[";
        for (auto taxon : taxa)
        {
            ss << taxon;
            if (taxon != taxa.back())
                ss << ",";
        }
        ss << "]\n";
        return ss.str();
    }
};

// Upstream result collection as map. Since tuples are not hashable, it is converted into a string before hashing.
struct TUpstreamKey
{
    using THash = std::string;
    TUpstreamKey(uint64_t taxid_, TKmerID fwd_, TKmerID rev_) : taxid(taxid_), fwd(fwd_), rev(rev_)
    {
    }

    // String represenation, usable as dictionary key or direct csv output.
    THash to_string()
    {
        std::stringstream s;
        s << taxid << "," << fwd << "," << rev;
        return s.str();
    }

    uint64_t taxid;
    TKmerID fwd;
    TKmerID rev;
};

template<typename uint_type>
std::string bits2str(uint_type i);

// Store enumerated length combinations of two kmers in two 64 bit unsigned integers.
// Enumerations follow lexicographical ordering. Since a kmerID may store up to
// 10 different kmer lengths, we have 100 possible kmer combinations. One bit for
// each kmer combination is reserved in the mask in big endian fashion, s.t. mask
// can be seen as the concatenation of 2x64 bits.
// 0: 0 with 0, i.e. length pattern l_min of kmerID1 combined with l_min of kmerID2
// 1: 0 with 1
// x: x/10 with x%10
template<typename TKmerID, typename TKmerLength>
struct TCombinePattern
{
private:
    // TODO: use union type for mask when doing SIMD vectorization
    std::bitset<100> data;  //? or union with SIMD 128

public:

    using TOffset = uint8_t;

    // return true if at least one combination bit is set.
    inline bool is_set()
    {
        return data.any();
    }
    // set a kmer combination by its lengths given the maximal length difference
    inline void set(uint64_t const prefix1, uint64_t const prefix2) noexcept
    {
        auto idx = __builtin_clzl(prefix1) * PREFIX_SIZE + __builtin_clzl(prefix2); // in [0:l_max^2[
        data.set(idx);
    }

    // unset bit if length combination doesn't pass a filter anymore
    // To be reset bit is expressed as offset w.r.t. PRIMER_MIN_LEN.
    inline void reset(TOffset const k_offset1, TOffset const k_offset2)
    {
        assert(k_offset1 <= PREFIX_SIZE && k_offset2 <= PREFIX_SIZE);
        data.reset(k_offset1 * PREFIX_SIZE + k_offset2);
    }

    // The number of combinations stored in data.
    uint64_t size()
    {
        //std::cout << "Enter size in cp: popcount_0 = " << __builtin_popcountll(data[0]) << ", popcount_1 = " << __builtin_popcountll(data[1]) << std::endl;
        return __builtin_popcountll(data[0]) + __builtin_popcountll(data[1]);
    }

    // Return all enumerated length combinations translated into kmer length offsets, i.e.
    // the true kmer length can be retrieved by adding PRIMER_MIN_LEN.
    void get_combinations(std::vector<std::pair<TOffset, TOffset>> & combinations)
    {
        combinations.clear();
        for (uint8_t i = 0; i < PREFIX_SIZE * PREFIX_SIZE; ++i)
        {
            if (data[i])
            {
                std::pair<TOffset, TOffset> pair{i / PREFIX_SIZE, (i % PREFIX_SIZE)};
                combinations.push_back(pair);
            }
        }
    }

    constexpr bool operator[](std::size_t pos) const
    {
        return data[pos];
    }

    // Wrapper for bitset::none(), returns true if no bit is set, else false.
    constexpr bool none() const noexcept
    {
        return data.none();
    }

    std::string to_string() const noexcept
    {
        return data.to_string();
    }
};

/*
 * A pair of encoded kmers (kmer IDs) given by their indices in the kmer ID
 * container associated to a reference bit vector. Since kmer IDs may encode up
 * to |[primer_max_length : primer_min_length]| kmers, the combine pattern stores
 * which kmers are combined. See encoding scheme here: TCombinePattern.
 */
template<typename TCombinePattern>
struct TPair
{
    TPair() = default;
    TPair(uint64_t reference_, uint64_t r_fwd_, uint64_t r_rev_, TCombinePattern cp_) :
        reference(reference_), r_fwd(r_fwd_), r_rev(r_rev_), cp(cp_) {}
    //std::tuple<uint64_t, uint64_t, TCombinePattern<TKmerID, TKmerLength>> TPair;
    ~TPair() = default;

    // Reference identifier (equivalent to position in corpus).
    uint64_t reference;
    // Rank of forward kmer.
    uint64_t r_fwd;
    // Rank of reverse kmer.
    uint64_t r_rev;
    // Length combinations of kmers as bit mask.
    TCombinePattern cp;
};

template<typename TPair>
using TPairList = std::vector<TPair>;

// Type for storing unique kmer combinations and their frequencies
// TODO: replace with TResult and vector<TResult>
// typedef std::pair<uint32_t, std::tuple<uint64_t, uint64_t, uint64_t, uint64_t> > TPairFreq;
// typedef std::vector<TPairFreq> TPairFreqList;

} // priset
