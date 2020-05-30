#pragma once

#include "simple_types.hpp"

// GenMap is built with SeqAn2 datatypes which need to be handled by PriSeT
namespace priset
{
    using dna = seqan::Dna5;

    // The Burrows-Wheeler Transfrom type in FM index compution.
    using TBWTLen = uint64_t;

    // The FM Index configuration type.
    using TFMIndexConfig = TGenMapFastFMIndexConfig<TBWTLen>;

    //
    using TString = seqan::String<seqan::Dna, seqan::Alloc<>>;

    //
    using TStringSet = seqan::StringSet<TString, seqan::Owner<seqan::ConcatDirect<SizeSpec_<priset::TSeqNo, TSeqPos>>>>;

    // The index type. TBiIndexConfig is defined in genmap/src/common.hpp.
    using TIndex = seqan::Index<TStringSet, TBiIndexConfig<TFMIndexConfig> >;

    // The DNA sequence type.
    using TSeq = seqan::String<priset::dna>;

    // The location type defined by sequence ID and position.
    using TLocation = seqan::Pair<priset::TSeqNo, priset::TSeqPos>;

    // The map of k-mer locations.
    using TLocations = std::map<TLocation, std::pair<std::vector<TLocation >, std::vector<TLocation>>>;

    // The handler type for the library sequences.
    // Library sequences are concatenated into one corpus.
    using TDirectoryInformation = typename seqan::StringSet<seqan::CharString, seqan::Owner<seqan::ConcatDirect<>>> ;

    // container for fasta header lines
    using TSequenceNames = typename seqan::StringSet<seqan::CharString, seqan::Owner<seqan::ConcatDirect<>>>;

    // container for fasta sequence lengths
    using TSequenceLengths = typename seqan::StringSet<uint32_t>;

}
