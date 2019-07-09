#pragma once

#include <algorithm>
#include <bitset>
#include <cmath>
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
void frequency_filter(priset::io_cfg_type const & io_cfg, primer_cfg_type const & primer_cfg, TLocations const & locations, TKmerLocations & kmer_locations, TKmerMap & kmer_map, TSeqNo  const min_occ) //, TDirectoryInformation const & directoryInformation, TSequenceNames & sequenceNames, TSequenceLengths & sequenceLengths)
{
    // = seqan::Pair<TSeqNo, TSeqPos>
    using TLocationKey = typename TLocations::key_type;
    // uniqueness indirectly preserved by (SeqNo, SeqPos) if list sorted lexicographically
    std::cout << "LOGGING: #kmers before pre frequency filtering = " << length(locations) << std::endl;

    if (!length(locations))
        return;
    //TSeqNo seqno;
    TSeqPos seqpos;
    // current location row
    std::vector<TLocationKey> row;
    // unique k-mers for retrieving sequences
    std::vector<seqan::Pair<TSeqNo, TSeqPos>> lookup();
    //TSeqPos seqpos_prev = 0;
    typename TLocations::const_iterator aux  = locations.end();
    --aux;

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
            //std::cout << "(" << seqan::getValueI1<TSeqNo, TSeqPos>(pair) << ", " << seqan::getValueI2<TSeqNo, TSeqPos>(pair) << ") ";
            row.push_back(pair);
        }
        std::cout << std::endl;
        // store locations,ID updated later
        kmer_locations.push_back(TKmerLocation{ID++, row}); //std::make_pair(ID++, row));
        row.clear();
    }
    lookup_sequences2(kmer_locations, kmer_map, io_cfg, primer_cfg);
    std::cout << "LOGGING: #kmers after frequency filtering: " << kmer_locations.size() << std::endl;
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
        kmer = kmer_map.at(kmer_locations[i].get_kmer_ID1());
        // filter by melting temperature
        if (kmer.Tm >= Tm_min && kmer.Tm <= Tm_max)
        {
            // filter by CG content
            if (filter_CG(primer_cfg, kmer.seq))
            {
                // Filter if Gibb's free energy is below -6 kcal/mol
                if (filter_self_dimerization(kmer.seq))
                    mask.set(i);
            }
        }
    }

    // Delete all masked out entries (mask_i = 0).
    for (int32_t i = kmer_locations.size() - 1; i >= 0; --i)
    {
        if (!mask[i])
        {
            kmer_map.erase(kmer_locations[i].get_kmer_ID()); // delete from dictionary
            kmer_locations.erase(kmer_locations.begin() + i); // erase associated locations
        }
    }
    std::cout << "LOGGING: # KMERS after chemical_filter_single: " << kmer_locations.size() << std::endl;
}

// check cross-dimerization.
void chemical_filter_pairs(/*primer_cfg_type const & primer_cfg, */TKmerPairs & kmer_pairs, TKmerMap & kmer_map)
{
    assert(kmer_pairs.size() < (1 << 12));
    std::bitset<1 << 12> mask{};
    uint16_t i = 0;
    for (auto kmer_pair : kmer_pairs)
    {
        if (filter_cross_dimerization(kmer_map[kmer_pair.get_kmer_ID1()], kmer_map[kmer_pair.get_kmer_ID2()]))
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
            unpaired.erase(kmer_pairs[i].get_kmer_ID1()), unpaired.erase(kmer_pairs[i].get_kmer_ID2());
        else
            kmer_pairs.erase(kmer_pairs.begin() + i);
    }
    // delete kmers from map that remained unpaired
    for (auto ID : unpaired)
        kmer_map.erase(ID);
    std::cout << "LOGGING: # KMER PAIRS after chemical filtering: " << kmer_pairs.size() << std::endl;
}

// post-filter candidates fulfilling chemical constraints by their relative frequency
/*void post_frequency_filter(TKmerLocations kmer_locations, TSeqNo occurrence_freq)
{

}*/

// filter k-mers by frequency and chemical properties
//template<typename TSequenceNames, typename TSequenceLengths>
void pre_filter_main(io_cfg_type const & io_cfg, primer_cfg_type const & primer_cfg, TLocations const & locations, TKmerLocations & kmer_locations, TKmerMap & kmer_map, TDirectoryInformation const & directoryInformation) //, TSequenceNames & sequenceNames, TSequenceLengths & sequenceLengths)
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
    // continue here
    min_occ = 1;
    std::cout << "MESSAGE: Cut-off frequency:\t" << min_occ << std::endl;
    // template<typename primer_cfg_type, typename TSeqSize, typename TSequenceNames, typename TSequenceLengths>
    // frequency filter and sequence fetching
    frequency_filter<TSequenceNames, TSequenceLengths>(io_cfg, primer_cfg, locations, kmer_locations, kmer_map, min_occ); //, directoryInformation, sequenceNames, sequenceLengths);
    print_kmer_locations(kmer_locations, kmer_map);
    chemical_filter_single(primer_cfg, kmer_locations, kmer_map);
//    chemical_filter_pairs();
    //post_frequency_filter(kmer_locations, primer_cfg.get_occurence_freq());

}

// combine helper, forward window until both iterators point to same reference ID or end
void fast_forward(TKmerLocation::TLocationVec const & locations1, TKmerLocations::value_type::const_iterator & it1_loc_start, TKmerLocation::TLocationVec const & locations2, TKmerLocations::value_type::const_iterator & it2_loc_start)
{
    std::cout << "enter fast_forward ..." << std::endl;
    while (it1_loc_start != locations1.end() && it2_loc_start != locations2.end() &&
            seqan::getValueI1<TSeqNo, TSeqPos>(*it1_loc_start) != seqan::getValueI1<TSeqNo, TSeqPos>(*it2_loc_start))
    {
        // speedup: use lower_bound/upper_bound
        while (it1_loc_start != locations1.end() && seqan::getValueI1<TSeqNo, TSeqPos>(*it1_loc_start) < seqan::getValueI1<TSeqNo, TSeqPos>(*it2_loc_start))
            ++it1_loc_start;
        while (it2_loc_start != locations2.end() && seqan::getValueI1<TSeqNo, TSeqPos>(*it2_loc_start) < seqan::getValueI1<TSeqNo, TSeqPos>(*it1_loc_start))
            ++it2_loc_start;
    }
    std::cout << "... leaving ff" << std::endl;
}

/* Combine based on suitable location distances s.t. transcript length is in permitted range.
 * Chemical suitability will be tested by a different function. First position indicates,
 * that the k-mer corresponds to a forward primer, and second position indicates reverse
 * primer, i.e. (k1, k2) != (k2, k1).
 */
void combine(primer_cfg_type const & primer_cfg, TKmerLocations const & kmer_locations, TKmerMap & kmer_map, TKmerPairs & kmer_pairs)
{
    //primer_cfg_type::size_interval_type transcript_range = primer_cfg.get_transcript_range();
    using it_loc_type = TKmerLocations::value_type::const_iterator;
    it_loc_type it1_loc_start, it1_loc_aux;
    it_loc_type it2_loc_start, it2_loc_aux;
    for (auto it1 = kmer_locations.begin(); it1 != kmer_locations.end() && it1 != kmer_locations.end()-1; ++it1)
    {
        for (auto it2 = it1+1; it2 != kmer_locations.end(); ++it2)
        {
            std::cout << "i1\n";
            // iterator to start position of current location for k-mer 1
            it1_loc_start = (*it1).locations.begin();
            // iterator to start position of current location for k-mer 2
            it2_loc_start = (*it2).locations.begin();
            // forward iterators to correspond to refer to same sequence ID or end
            fast_forward((*it1).locations, it1_loc_start, (*it2).locations, it2_loc_start);
            std::cout << "i2\n";
            // no common reference sequences => forward kmer iterators
            if (it1_loc_start == (*it1).locations.end() || it2_loc_start == (*it2).locations.end())
                continue;

            // foward iterator for k-mer 1 on same sequence (loc)
            it1_loc_aux = it1_loc_start;
            // foward iterator for k-mer 2 on same sequence (loc)
            it2_loc_aux = it2_loc_start;
            // invariant after entering this loop: loc(it1) == loc(it2)
            auto seq_ID = seqan::getValueI1<TSeqNo, TSeqPos>(*it1_loc_start);
            //for (it1_loc_aux = it1_loc_start, it2_loc_aux = it2_loc_start; )
            while (it1_loc_start != (*it1).locations.end() && it1_loc_aux != (*it1).locations.end() && it2_loc_start != (*it2).locations.end() && seq_ID == seqan::getValueI1<TSeqNo, TSeqPos>(*it2_loc_aux))
            {
//                assert((seqan::getValueI1<TSeqNo, TSeqPos>(*it1_loc_aux) == seqan::getValueI1<TSeqNo, TSeqPos>(*it2_loc_aux)));
                // valid combination?
                std::cout << "h4\n";
                auto pos_kmer1 = seqan::getValueI2<TSeqNo, TSeqPos>(*it1_loc_aux);
                if (it2_loc_aux != (*it2).locations.end())
                {
                    auto pos_kmer2 = seqan::getValueI2<TSeqNo, TSeqPos>(*it2_loc_aux);
                    auto pos_delta = (pos_kmer1 < pos_kmer2) ? pos_kmer2 - pos_kmer1 : pos_kmer1 - pos_kmer2;
                    pos_delta -= primer_cfg.get_kmer_length();
                    if (pos_delta >= primer_cfg.get_transcript_range().first && pos_delta <= primer_cfg.get_transcript_range().second)
                    {
                        TKmerID kmer_fwd_new = (pos_kmer1 < pos_kmer2) ? (*it1).get_kmer_ID() : (*it2).get_kmer_ID();
                        TKmerID kmer_rev_new = (pos_kmer1 < pos_kmer2) ? (*it2).get_kmer_ID() : (*it1).get_kmer_ID();
                        auto pair_location = std::make_tuple(seq_ID, std::min<TSeqPos>(pos_kmer1, pos_kmer2), std::max<TSeqPos>(pos_kmer1, pos_kmer2));
                        std::cout << "(kmer_fwd_new, kmer_rev_new) = (" << kmer_fwd_new << ", " << kmer_rev_new << ") at [(refID = " << std::get<0>(pair_location) << " at positions: " << std::get<1>(pair_location) << ", " << std::get<2>(pair_location) << ")]\n";
                        // extend location vector if pair combinations already in result
                        if (kmer_pairs.size() && kmer_pairs.back().get_kmer_ID1() == kmer_fwd_new && kmer_pairs.back().get_kmer_ID2() == kmer_rev_new)
                        {
                            std::cout << "j1\n";
                            kmer_pairs[kmer_pairs.size()-1].pair_locations.push_back(pair_location);
                        }
                        else
                        {
                            std::cout << "j1\n";
                            //TKmerPair::TKmerPairLocations first_pair{pair_location};
                            TKmerPair pair{kmer_fwd_new, kmer_rev_new, abs(kmer_map.at(kmer_fwd_new).Tm - kmer_map.at(kmer_rev_new).Tm), pair_location};
                            kmer_pairs.push_back(pair);
                        }
                        std::cout << "h1\n";
                    }
                    ++it2_loc_aux;
                }
                std::cout << "h2\n";
                // all combinations tested for second k-mer
                if (it2_loc_aux == (*it2).locations.end())
                {
                    std::cout << "LOGGING: forward it1 and reset it2\n";
                    // reset to start of current sequence if next k-mer of it1 refers to same sequence, else forward it2
                    if (++it1_loc_aux == (*it1).locations.end())
                        break;
                    if (seqan::getValueI1<TSeqNo, TSeqPos>(*it1_loc_aux) == pos_kmer1)
                        it2_loc_aux = it2_loc_start;
                    else
                    {
                        it1_loc_start = it1_loc_aux;
                        fast_forward((*it1).locations, it1_loc_start, (*it2).locations, it2_loc_start);
                        it1_loc_aux = it1_loc_start;
                        it2_loc_aux = it2_loc_start;
                    }
                }
                std::cout << "h3\n";
            }
        }
    }
    std::cout << "LOGGING: # KMER PAIRS combined from " << kmer_locations.size() << " kmers: " << kmer_pairs.size() << std::endl;
}

/*void post_filter_main(primer_cfg_type const & primer_cfg, TKmerLocations & kmer_locations, TKmerPairs & pairs)
{

}*/

}  // namespace priset
