#pragma once

#include <algorithm>
#include <bitset>
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
template<typename TSeqSize, typename TSequenceNames, typename TSequenceLengths>
void frequency_filter(priset::io_cfg_type const & io_cfg, primer_cfg_type const & primer_cfg, TLocations const & locations, TKmerLocations & kmer_locations, TSeqSize  const min_occ, TDirectoryInformation const & directoryInformation, TSequenceNames & sequenceNames, TSequenceLengths & sequenceLengths)
{
    // = seqan::Pair<TSeqNo, TSeqPos>
    using TLocationKey = typename TLocations::key_type;
    // uniqueness indirectly preserved by (SeqNo, SeqPos) if list sorted lexicographically
    std::cout << "LOGGING: number of kmers before pre frequency filtering = " << length(locations) << std::endl;

    std::cout << "unique k-mers:\n";
    //TSeqNo seqno;
    TSeqPos seqpos;
    // current location row
    std::vector<TLocationKey> row;
    // unique k-mers for retrieving sequences
    std::vector<seqan::Pair<TSeqNo, TSeqPos>> lookup();
    for (typename TLocations::const_iterator it = locations.begin(); it != locations.end(); ++it)
    {
        //seqno = seqan::getValueI1<TSeqNo, TSeqPos>(it->first);
        seqpos = seqan::getValueI2<TSeqNo, TSeqPos>(it->first);
        // not enough k-mer occurences => continue
        if ((it->second).first.size() < min_occ)
            continue;
        // use symmetry and lexicographical ordering of locations to skip already seen ones
        if (seqan::getValueI1<TSeqNo, TSeqPos>((it->second).first[0]) < seqpos)
            continue;
        // invariant: min_occ is always ≥ 2
        for (seqan::Pair<TSeqNo, TSeqPos> pair : (it->second).first)
        {
            std::cout << "(" << seqan::getValueI1<TSeqNo, TSeqPos>(pair) << ", " << seqan::getValueI2<TSeqNo, TSeqPos>(pair) << ") ";
            row.push_back(pair);
        }
        std::cout << std::endl;
        // store locations, dna sequence retrieved later
        kmer_locations.push_back(std::make_pair(TSeq(0), row));
        row.clear();
    }
    lookup_sequences<primer_cfg_type>(kmer_locations, io_cfg, primer_cfg, directoryInformation);
    std::cout << "LOGGING: number of kmers before after frequency filtering = " << kmer_locations.size() << std::endl;
}

/*
 * Filter kmers based on their chemical properties regardless of their pairing.
 * Constraints that are checked: melting tempaerature, CG content
 */
 void chemical_filter_single(primer_cfg_type const & primer_cfg, TKmerLocations & kmer_locations)
 {
     std::cout << "LOGGING: number of kmers before chemical_filter_single = " << kmer_locations.size() << std::endl;
     //typedef std::vector<std::pair<TSeq, std::vector<seqan::Pair<priset::TSeqNo, priset::TSeqPos> > > > TKmerLocations;
     assert(kmer_locations.size() < (1 << 18));
     seqan::String<priset::dna> kmer_seq;
     auto mask = std::bitset<1 << 18>{}.set();

     // Filter by melting temperature
     std::bitset<1 << 18> aux_filter{};
     for (TKmerLocations::size_type i = 0; i < kmer_locations.size(); ++i)
     {
         kmer_seq = kmer_locations[i].first;
         if (chemistry::filter_Tm(primer_cfg, kmer_seq))
            aux_filter.set(i);
     }
    mask &= aux_filter, aux_filter.reset();

    // Filter by CG content
    for (TKmerLocations::size_type i = 0; i < kmer_locations.size(); ++i)
    {
        if (mask[i])
        {
            kmer_seq = kmer_locations[i].first;
            if (chemistry::filter_CG<primer_cfg_type>(primer_cfg, kmer_seq))
                aux_filter.set(i);
        }
    }
    mask &= aux_filter, aux_filter.reset();

    // Filter if Gibb's free energy is below -6 kcal/mol
    for (TKmerLocations::size_type i = 0; i < kmer_locations.size(); ++i)
    {
        if (mask[i])
        {
            kmer_seq = kmer_locations[i].first;
            if (chemistry::filter_self_dimerization(kmer_seq))
                aux_filter.set(i);
        }
    }
    mask &= aux_filter;

    // Delete all masked out entries (mask_i = 0). Vector erasure invalidates iterators
    // and references at or after the point of the erase, including the end() iterator,
    // that's why we erase from back to front.
    for (auto i = kmer_locations.size()-1; i >= 0; --i)
    {
        if (!mask[i])
            kmer_locations.erase(kmer_locations.begin() + i, kmer_locations.begin() + i + 1);
    }

     std::cout << "LOGGING: number of kmers after chemical_filter_single = " << kmer_locations.size() << std::endl;
}

// check cross-dimerization.
void chemical_filter_pairs(primer_cfg_type const & primer_cfg, TKmerLocations & kmer_locations, TMatrix & kmer_pairs)
{
    for (uint64_t i = 0; i < )
}

// post-filter candidates fulfilling chemical constraints by their relative frequency
void post_frequency_filter(TKmerLocations kmer_locations, TSeqNo occurrence_freq)
{

}

// filter k-mers by frequency and chemical properties
template<typename TSequenceNames, typename TSequenceLengths>
void pre_filter_main(io_cfg_type const & io_cfg, primer_cfg_type const & primer_cfg, TLocations const & locations, TKmerLocations & kmer_locations, TDirectoryInformation const & directoryInformation, TSequenceNames & sequenceNames, TSequenceLengths & sequenceLengths)
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
    // template<typename primer_cfg_type, typename TSeqSize, typename TSequenceNames, typename TSequenceLengths>
    frequency_filter<primer_cfg_type, TSeqNo, TSequenceNames, TSequenceLengths>(io_cfg, primer_cfg, locations, kmer_locations, min_occ, directoryInformation, sequenceNames, sequenceLengths);
    chemical_filter_single<primer_cfg_type>(primer_cfg, kmer_locations);
//    chemical_filter_pairs();
    //post_frequency_filter(kmer_locations, primer_cfg.get_occurence_freq());

}


/* Combine based on suitable location distances s.t. transcript length is in permitted range.
 * Chemical suitability will be tested.
 * Details: for each combination of two different kmers which are associated with sorted
 * location lists (firstly by ID, secondly by sequence position).
 * TODO: verify sorting by ID (1) and sequence position (2)
 */
void combine(primer_cfg_type const & primer_cfg, TKmerLocations & kmer_locations, TPairs & pairs)
{
    auto transcript_range = primer_cfg.get_transcript_range();
    for (auto it1 = kmer_locations.begin(); it1 != kmer_locations.end()-1; ++it1)
    {
        for (auto it2 = it1+1; it2 != kmer_locations.end(); ++it2)
        {
            // parse location lists linearly for each kmer combination
            auto it1_loc = (*it1).second.begin();
            auto it2_loc = (*it2).second.begin();
            while (it1_loc != (*it1).second.end() && it2_loc != (*it2).second.end())
            {
                // Continue here
            }
        }
    }
    pairs.resize();
}

// Filter pairs
void post_filter_main(TKmerLocations & kmer_locations, TMatrix & pairs)
{

}

}  // namespace priset
