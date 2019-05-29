#pragma once

#include <algorithm>
#include <bitset>
#include <fstream>
#include <string>
#include <unordered_set>

#include "../submodules/genmap/src/common.hpp"
#include "../submodules/genmap/src/genmap_helper.hpp"

#include "primer_cfg_type.hpp"
#include "types.hpp"
#include "utilities.hpp"

namespace priset
{
// TODO: globally or hierarchical?
// pre-filter and sequence fetch
// 1. filter candidates by number of occurences only independent of their chemical suitability
// 2. fetch sequence and check chemical constraints that need to hold for a single primer
// frequency_filter<primer_cfg_type, TSeqNo, TSequenceNames, TSequenceLengths>(io_cfg, primer_cfg, locations, kmer_locations, min_occ, directoryInformation, sequenceNames, sequenceLengths);
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
    TSeqPos seqpos_prev = 0;
    for (typename TLocations::const_iterator it = locations.begin(); it != locations.end(); ++it)
    {
        //seqno = seqan::getValueI1<TSeqNo, TSeqPos>(it->first);
        seqpos = seqan::getValueI2<TSeqNo, TSeqPos>(it->first);
        if (seqpos < seqpos_prev)
        {
            std::cout << "seqpos < seqpos_prev: " << seqpos << ", " << seqpos_prev << std::endl;
            exit(0);
        }
        seqpos_prev = seqpos;
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
        kmer_locations.push_back(std::make_pair(TKmer{}, row));
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
    assert(kmer_locations.size() < (1 << 18));
    auto mask = std::bitset<1 << 18>{}.set();
    std::bitset<1 << 18> aux_filter{};
    // Filter by melting temperature
    float Tm_min = primer_cfg.get_min_Tm();
    float Tm_max = primer_cfg.get_max_Tm();
    TKmer kmer;
    for (TKmerLocations::size_type i = 0; i < kmer_locations.size(); ++i)
    {
        kmer = kmer_locations[i].first;
        if (kmer.Tm >= Tm_min && kmer.Tm <= Tm_max)
            aux_filter.set(i);
    }
    mask &= aux_filter, aux_filter.reset();

    // Filter by CG content
    for (TKmerLocations::size_type i = 0; i < kmer_locations.size(); ++i)
    {
        if (mask[i])
        {
            if (chemistry::filter_CG<primer_cfg_type>(primer_cfg, kmer_locations[i].first.seq))
                aux_filter.set(i);
        }
    }
    mask &= aux_filter, aux_filter.reset();

    // Filter if Gibb's free energy is below -6 kcal/mol
    for (TKmerLocations::size_type i = 0; i < kmer_locations.size(); ++i)
    {
        if (mask[i])
        {
            if (chemistry::filter_self_dimerization(kmer_locations[i].first.seq))
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
void chemical_filter_pairs(primer_cfg_type const & primer_cfg, TKmerLocations & kmer_locations, TPairs & kmer_pairs)
{
    assert(kmer_pairs.size() < (1 << 12));
    std::bitset<1 << 12> mask{};
    uint16_t i = 0;
    for (auto pair : kmer_pairs)
    {
        if (chemistry::filter_cross_dimerization(kmer_locations[pair.first].first.seq, kmer_locations[pair.second].first.seq))
            mask.set(i);
        ++i;
    }

    std::unordered_set<TKmerID> paired{};
    // delete unset entries
    for (auto i = kmer_pairs.size()-1; i >= 0; --i)
    {
        if (mask[i])
            paired.insert(kmer_pairs[i].first), paired.insert(kmer_pairs[i].second);
        else
            kmer_pairs.erase(kmer_pairs.begin() + i);
    }
    // delete kmers that are unpaired
    for (uint16_t i = kmer_locations.size(); i >= 0; --i)
    {
        if (paired.find(kmer_locations[i].first.ID) == paired.end())
            kmer_locations.erase(kmer_locations.begin() + i);
    }
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
    frequency_filter<TSeqNo, TSequenceNames, TSequenceLengths>(io_cfg, primer_cfg, locations, kmer_locations, min_occ, directoryInformation, sequenceNames, sequenceLengths);
    chemical_filter_single(primer_cfg, kmer_locations);
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
    primer_cfg_type::size_interval_type transcript_range = primer_cfg.get_transcript_range();
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
}

// Filter pairs

/*void post_filter_main(primer_cfg_type const & primer_cfg, TKmerLocations & kmer_locations, TPairs & pairs)
{

}*/

}  // namespace priset
