// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================
//          Author: Marie Hoffmann <ozymandiaz147 AT gmail.com>
//          Manual: https://github.com/mariehoffmann/PriSeT

#pragma once

// #include <stdlib.h>
#include <vector>

// #include "utilities.hpp"

#include "../../submodules/sdsl-lite/include/sdsl/bit_vectors.hpp"
#include "../../submodules/genmap/src/common.hpp"

#define ONE_LSHIFT_63 (1ULL << 63)

namespace priset
{

enum TIMEIT {
    MAP, // mappability computation
    FILTER1_TRANSFORM, //
    COMBINE_FILTER2, // kmer combiner
    PAIR_FREQ, // Pair frequency cutoff
    SIZE
};

// type declarations
using TSeqNo = uint64_t;
using TSeqPos = uint64_t;
// Kmer length type. A negative indicates reverse direction given associated position.
using TKmerLength = int64_t;

// The type for storing up to 10 k-mers which are prefixes from each other.
typedef uint64_t TKmerID;

// static TKmerID NULL_TKMERID = 0;

// The taxonomic identifier type.
using Taxid = uint64_t;

// The type of an accession identifier, e.g. from NCBI GenBank.
using Accession = std::string;

// The type for identifying numerically an accession.
using AccessionID = uint32_t;

// The type for numerical accession identifiers, 1-based.
// typedef uint64_t AccessionID;

// The type of an accession
// typedef std::string Accession;

// A bit vector for each reference, with 1 indicating a kmer starting position.
typedef sdsl::bit_vector TReference;
typedef std::vector<TReference> TReferences;

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

} // priset
