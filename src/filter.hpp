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
template<typename TSequenceNames, typename TSequenceLengths>
void frequency_filter(priset::io_cfg_type const & io_cfg, primer_cfg_type const & primer_cfg, TLocations const & locations, TKmerLocations & kmer_locations, TKmerMap & kmer_map, TSeqNo  const min_occ, TDirectoryInformation const & directoryInformation, TSequenceNames & sequenceNames, TSequenceLengths & sequenceLengths)
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
    typename TLocations::const_iterator aux  = locations.end();
    --aux;
    std::cout << "last unfiltered location: " << seqan::getValueI2<TSeqNo, TSeqPos>(aux->first);
    for (seqan::Pair<TSeqNo, TSeqPos> pair : aux->second.first)
    {
        std::cout << "(" << seqan::getValueI1<TSeqNo, TSeqPos>(pair) << ", " << seqan::getValueI2<TSeqNo, TSeqPos>(pair) << ") ";
    }
    //exit(0);
    TKmerID ID = 1;
    for (typename TLocations::const_iterator it = locations.begin(); it != locations.end(); ++it)
    {
        seqpos = seqan::getValueI2<TSeqNo, TSeqPos>(it->first);
        // not enough k-mer occurences => continue
        if ((it->second).first.size() < min_occ)
            continue;
        // use symmetry and lexicographical ordering of locations to skip already seen ones
        if (it->second.first.size() && seqan::getValueI1<TSeqNo, TSeqPos>(it->second.first[0]) < seqpos)
            continue;
        // invariant: min_occ is always ≥ 2
        for (seqan::Pair<TSeqNo, TSeqPos> pair : it->second.first)
        {
            std::cout << "(" << seqan::getValueI1<TSeqNo, TSeqPos>(pair) << ", " << seqan::getValueI2<TSeqNo, TSeqPos>(pair) << ") ";
            row.push_back(pair);
        }
        std::cout << std::endl;
        // store locations,ID updated later
        kmer_locations.push_back(std::make_pair(ID++, row));
        row.clear();
    }
    lookup_sequences<primer_cfg_type>(kmer_locations, kmer_map, io_cfg, primer_cfg, directoryInformation);
    std::cout << "LOGGING: # KMERS after frequency filtering: " << kmer_locations.size() << std::endl;
}

/*
 * Filter kmers based on their chemical properties regardless of their pairing.
 * Constraints that are checked: melting tempaerature, CG content
 */
void chemical_filter_single(primer_cfg_type const & primer_cfg, TKmerLocations & kmer_locations, TKmerMap & kmer_map)
{
    assert(kmer_locations.size() < (1 << 18));
    auto mask = std::bitset<1 << 18>{};
    //std::bitset<1 << 18> aux_filter{};
    // Filter by melting temperature
    float Tm_min = primer_cfg.get_min_Tm();
    float Tm_max = primer_cfg.get_max_Tm();
    TKmer kmer;
    for (TKmerLocations::size_type i = 0; i < kmer_locations.size(); ++i)
    {
        kmer = kmer_map.at(kmer_locations[i].first);
        // filter by melting temperature
        if (kmer.Tm >= Tm_min && kmer.Tm <= Tm_max)
        {
            // filter by CG content
            if (chemistry::filter_CG<primer_cfg_type>(primer_cfg, kmer.seq))
            {
                // Filter if Gibb's free energy is below -6 kcal/mol
                if (chemistry::filter_self_dimerization(kmer.seq))
                    mask.set(i);
            }
        }
    }

    // Delete all masked out entries (mask_i = 0).
    for (int32_t i = kmer_locations.size() - 1; i >= 0; --i)
    {
        if (!mask[i])
        {
            kmer_map.erase(kmer_locations[i].first); // delete from dictionary
            kmer_locations.erase(kmer_locations.begin() + i); // erase associated locations
        }
    }

     std::cout << "LOGGING: # KMERS after chemical_filter_single: " << kmer_locations.size() << std::endl;
}

// check cross-dimerization.
void chemical_filter_pairs(primer_cfg_type const & primer_cfg, TKmerPairs & kmer_pairs, TKmerMap & kmer_map)
{
    assert(kmer_pairs.size() < (1 << 12));
    std::bitset<1 << 12> mask{};
    uint16_t i = 0;
    for (auto kmer_pair : kmer_pairs)
    {
        if (chemistry::filter_cross_dimerization(kmer_map[kmer_pair.kmer_fwd], kmer_map[kmer_pair.kmer_rev]))
            mask.set(i);
        ++i;
    }

    std::unordered_set<TKmerID> unpaired{};
    for (auto it = kmer_map.begin(); it != kmer_map.end(); ++it)
        unpaired.insert(it->first);
    // delete unset entries
    for (int16_t i = kmer_pairs.size()-1; i >= 0; --i)
    {
        if (mask[i])
            unpaired.erase(kmer_pairs[i].kmer_fwd), unpaired.erase(kmer_pairs[i].kmer_rev);
        else
            kmer_pairs.erase(kmer_pairs.begin() + i);
    }
    // delete kmers from map that remained unpaired
    for (auto ID : unpaired)
        kmer_map.erase(ID);
    std::cout << "LOGGING: # KMER PAIRS after chemical filtering: " << kmer_pairs.size() << std::endl;
}

// post-filter candidates fulfilling chemical constraints by their relative frequency
void post_frequency_filter(TKmerLocations kmer_locations, TSeqNo occurrence_freq)
{

}

// filter k-mers by frequency and chemical properties
template<typename TSequenceNames, typename TSequenceLengths>
void pre_filter_main(io_cfg_type const & io_cfg, primer_cfg_type const & primer_cfg, TLocations const & locations, TKmerLocations & kmer_locations, TKmerMap & kmer_map, TDirectoryInformation const & directoryInformation, TSequenceNames & sequenceNames, TSequenceLengths & sequenceLengths)
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
    // frequency filter and sequence fetching
    frequency_filter<TSequenceNames, TSequenceLengths>(io_cfg, primer_cfg, locations, kmer_locations, kmer_map, min_occ, directoryInformation, sequenceNames, sequenceLengths);
    print_kmer_locations(kmer_locations, kmer_map);
    chemical_filter_single(primer_cfg, kmer_locations, kmer_map);
//    chemical_filter_pairs();
    //post_frequency_filter(kmer_locations, primer_cfg.get_occurence_freq());

}


/* Combine based on suitable location distances s.t. transcript length is in permitted range.
 * Chemical suitability will be tested by a different function. First position indicates,
 * that the k-mer corresponds to a forward primer, and second position indicates reverse
 * primer, i.e. (k1, k2) != (k2, k1).
 */
void combine(primer_cfg_type const & primer_cfg, TKmerLocations & kmer_locations, TKmerPairs & pairs)
{
    primer_cfg_type::size_interval_type transcript_range = primer_cfg.get_transcript_range();
    using it_loc_type = TKmerLocations::value_type::second_type::const_iterator;
    it_loc_type it1_start, it1_aux, it2_start, it2_aux;
    for (auto it1 = kmer_locations.begin(); it1 != kmer_locations.end()-1; ++it1)
    {
        for (auto it2 = it1+1; it2 != kmer_locations.end(); ++it2)
        {
            // iterator to start position of current location for k-mer 1
            it1_start = (*it1).second.begin();
            // iterator to start position of current location for k-mer 2
            it2_start = (*it2).second.begin();
            // forward iterators to correspond to refer to same sequence ID or end
            while (seqan::getValueI1<TSeqNo, TSeqPos>(*it1_start) != seqan::getValueI1<TSeqNo, TSeqPos>(*it2_start))
            {
                while (it1_start != (*it1).second.end() && seqan::getValueI1<TSeqNo, TSeqPos>(*it1_start) < seqan::getValueI1<TSeqNo, TSeqPos>(*it2_start))
                    ++it1_start;
                while (it2_start != (*it2).second.end() && seqan::getValueI1<TSeqNo, TSeqPos>(*it2_start) < seqan::getValueI1<TSeqNo, TSeqPos>(*it1_start))
                    ++it2_start;
            }
            // foward iterator for k-mer 1 on same sequence (loc)
            it1_aux = it1_start;
            // foward iterator for k-mer 2 on same sequence (loc)
            it2_aux = it2_start;
            // invariant after entering this loop: loc(it1) == loc(it2)
            while (it1_start != (*it1).second.end() && it1_aux != (*it1).second.end() && it2_start != (*it2).second.end())
            {
                // valid combination?
                auto pos_kmer1 = seqan::getValueI2<TSeqNo, TSeqPos>(*it1_aux);
                auto pos_kmer2 = seqan::getValueI2<TSeqNo, TSeqPos>(*it2_aux);
                if (it2_aux != (*it2).second.end())
                {
                    auto pos_delta = (pos_kmer1 < pos_kmer2) ? pos_kmer2 - pos_kmer1 : pos_kmer1 - pos_kmer2;
                        pos_delta += primer_cfg.get_kmer_length();
                        if (pos_delta >= primer_cfg.get_transcript_range().first && pos_delta <= primer_cfg.get_transcript_range().second)
                            pairs.push_back((pos_kmer1 < pos_kmer2) ? TPair{(*it1).first, (*it2).first} : TPair{(*it2).first, (*it1).first});
                }
                // all combinations tested for second k-mer
                if (it2_aux == (*it2).second.end())
                {
                    // reset to start of current sequence if next k-mer of it1 refers to same sequence, else forward it2
                    if (++it1_start == (*it1).second.end())
                        break;
                    while (seqan::getValueI1<TSeqNo, TSeqPos>(*it1_start) != seqan::getValueI1<TSeqNo, TSeqPos>(*++it2_start));
                    it2_aux = it2_start;

                }
            }
        }
    }
    std::cout << "LOGGING: # KMER PAIRS combined from " << kmer_locations.size() << " kmers: " << pairs.size() << std::endl;
}

/*void post_filter_main(primer_cfg_type const & primer_cfg, TKmerLocations & kmer_locations, TKmerPairs & pairs)
{

}*/

}  // namespace priset
