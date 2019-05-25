#pragma once

#include <algorithm>
#include <fstream>
#include <string>

#include "../submodules/genmap/src/common.hpp"
#include "../submodules/genmap/src/genmap_helper.hpp"

#include "types.hpp"
#include "utilities.hpp"

namespace priset
{
// TODO: globally or hierarchical?
// pre-filter and sequence fetch
// 1. filter candidates by number of occurences only independent of their chemical suitability
// 2. fetch sequence and check chemical constraints that need to hold for a single primer
template<typename primer_cfg_type, typename TLocations, typename TKmerLocations, typename TSeqSize, typename TDirectoryInformation, typename TSequenceNames, typename TSequenceLengths>
void pre_frequency_filter(priset::io_config & io_cfg, primer_cfg_type & primer_cfg, TLocations & locations,
    TKmerLocations & kmer_locations, TSeqSize min_occ, TDirectoryInformation & directoryInformation, TSequenceNames & sequenceNames, TSequenceLengths & sequenceLengths)
{
    // = seqan::Pair<TSeqNo, TSeqPos>
    using TLocationKey = typename TLocations::key_type;
    // uniqueness indirectly preserved by (SeqNo, SeqPos) if list sorted lexicographically

    std::cout << "unique k-mers:\n";
    TSeqNo seqno;
    TSeqPos seqpos;
    // current location row
    std::vector<TLocationKey> row;
    // unique k-mers for retrieving sequences
    std::vector<seqan::Pair<TSeqNo, TSeqPos>> lookup();
    for (typename TLocations::const_iterator it = locations.begin(); it != locations.end(); ++it)
    {
        seqno = seqan::getValueI1<TSeqNo, TSeqPos>(it->first);
        seqpos = seqan::getValueI2<TSeqNo, TSeqPos>(it->first);
        // not enough k-mer occurences => continue
        if ((it->second).first.size() < min_occ)
            continue;
        // use symmetry and lexicographical ordering of locations to skip already seen ones
        if (seqan::getValueI1<TSeqNo, TSeqPos>((it->second).first[0]) < seqpos)
            continue
        // invariant: min_occ is always ≥ 2
        for (seqan::Pair<TSeqNo, TSeqPos> pair : (it->second).first)
        {
            std::cout << "(" << seqan::getValueI1<TSeqNo, TSeqPos>(pair) << ", " << seqan::getValueI2<TSeqNo, TSeqPos>(pair) << ") ";
            row.push_back(pair);
        }
        // store locations, dna sequence retrieved later
        kmer_locations.push_back(std::make_pair(TSeq(0), row));
        row.clear();
    }

    //lookup_sequences<primer_cfg_type>(kmer_locations, io_cfg, primer_cfg);
}

// post-filter candidates fulfilling chemical constraints by their relative frequency
template<typename TKmerLocations, typename TSeqNo>
void post_frequency_filter(TKmerLocations kmer_locations, TSeqNo occurrence_freq)
{

}

// filter k-mers by frequency and chemical properties
template<typename primer_cfg_type, typename TLocations, typename TKmerLocations, typename TDirectoryInformation, typename TSequenceNames, typename TSequenceLengths>
void filter(priset::io_config & io_cfg, primer_cfg_type & primer_cfg, TLocations & locations, TKmerLocations & kmer_locations, TDirectoryInformation & directoryInformation, TSequenceNames & sequenceNames, TSequenceLengths & sequenceLengths)
{
    using TSeqNo = typename seqan::Value<typename TLocations::key_type, 1>::Type;
    // get number of taxids with at least one accession assigned to it
    std::string acc_file_name = io_cfg.get_acc_file();
    std::string line;
    TSeqNo min_occ = 0;
    std::ifstream in;
    in.open(acc_file_name);
    while(!in.eof()) {
	       getline(in, line);
	       ++min_occ;
    }
    in.close();
    min_occ = (min_occ) ? min_occ - 1 : 0;
    std::cout << "STATS: Number of taxids with one or more accessions:\t" << min_occ << std::endl;
    // scale to be lower frequency bound for filters
    min_occ = std::max<TSeqNo>(2, TSeqNo(float(min_occ) * primer_cfg.get_occurence_freq()));
    std::cout << "MESSAGE: Cut-off frequency:\t" << min_occ << std::endl;
    pre_frequency_filter<primer_cfg_type, TLocations, TKmerLocations, TSeqNo, TDirectoryInformation, TSequenceNames, TSequenceLengths>(io_cfg, primer_cfg, locations, kmer_locations, min_occ, directoryInformation, sequenceNames, sequenceLengths);

    //post_frequency_filter(kmer_locations, primer_cfg.get_occurence_freq());

}

}  // namespace priset
