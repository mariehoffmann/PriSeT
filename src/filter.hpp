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
// TODO: kmer_location is the transformed to reference bit vector
// pre-filter and sequence fetch
// 1. filter candidates by number of occurences only independent of their chemical suitability
// 2. fetch sequence and check chemical constraints that need to hold for a single primer
void frequency_filter(priset::io_cfg_type const & io_cfg, TKLocations & locations, TKmerLocations & kmer_locations, TSeqNo const cutoff)
{
    // uniqueness indirectly preserved by (SeqNo, SeqPos) if list sorted lexicographically
    assert(length(locations));

    // load corpus
    seqan::StringSet<seqan::DnaString, seqan::Owner<seqan::ConcatDirect<>>> text;
    fs::path text_path = io_cfg.get_index_txt_path();
    std::cout << "text_path: " << text_path << std::endl;
    seqan::open(text, text_path.string().c_str(), seqan::OPEN_RDONLY);

    TKmerID kmer_ID;
    std::vector<TLocation> fwd;
    for (typename TKLocations::const_iterator it = locations.begin(); it != locations.end(); ++it)
    {
        // not enough k-mer occurences => continue
        // Note: getOccurrences and resetLimits in genmap lead to less kmers occurences than countOccurrences
        if ((it->second).first.size() < cutoff)
            continue;

        const auto & [seqNo, seqPos, K] = (it->first);

        // use symmetry and lexicographical ordering of locations to skip already seen ones
        if (it->second.first.size() && seqan::getValueI1<TSeqNo, TSeqPos>(it->second.first[0]) < seqPos)
            continue;
        // invariant: cutoff is always ≥ 2
        for (TLocation pair : it->second.first)
        {
            fwd.push_back(pair);
        }
        // encode dna sequence
        seqan::DnaString seq = seqan::valueById(text, seqNo);
        auto const & kmer_str = seqan::infixWithLength(seq, seqPos, K);

        kmer_ID = dna_encoder(kmer_str);
        // replace kmer_ID
        kmer_locations.push_back(TKmerLocation{kmer_ID, K, fwd});
        // TODO: same for reverse?
        fwd.clear();
    }
    locations.clear();
}

/*
 * Filter kmers based on their chemical properties regardless of their pairing.
 * Constraints that are checked: melting tempaerature, CG content
 */
void chemical_filter_single(primer_cfg_type const & primer_cfg, TKmerLocations & kmer_locations)
{
    assert(kmer_locations.size() < (1 << 24));
    auto mask = std::bitset<1 << 24>{};
    //std::bitset<1 << 18> aux_filter{};
    // Filter by melting temperature
    float Tm_min = primer_cfg.get_min_Tm();
    float Tm_max = primer_cfg.get_max_Tm();
    TKmerID kmer_ID;

    for (TKmerLocations::size_type i = 0; i < kmer_locations.size(); ++i)
    {

        // TODO: optimize - drop kmer when computing Tm in sequence lookup fct
        kmer_ID = kmer_locations[i].get_kmer_ID1();
        //std::cout << "kmer_ID = " << kmer_ID << std::endl;
        auto Tm = primer_melt_wallace(kmer_ID);
    //    std::cout << "Tm = " << Tm << std::endl;
        // filter by melting temperature
        if (Tm >= Tm_min && Tm <= Tm_max)
        {
            // filter by CG content
    //        std::cout << "call filter_CG ...\n";
            if (filter_CG(primer_cfg, kmer_ID))
            {

                // Filter if Gibb's free energy is below -6 kcal/mol
                if (filter_self_dimerization(kmer_ID))
                {
                    if (filter_repeats_runs(kmer_ID))
                        mask.set(i);
                }
            }
        }
    }

    // Delete all masked out entries (mask_i = 0).
    for (int32_t i = kmer_locations.size() - 1; i >= 0; --i)
    {
        if (!mask[i])
            kmer_locations.erase(kmer_locations.begin() + i);
    }
}

// check cross-dimerization.
void chemical_filter_pairs(/*primer_cfg_type const & primer_cfg, */TKmerPairs & kmer_pairs)
{
    assert(kmer_pairs.size() < (1 << 24));
    std::bitset<1 << 24> mask{};
    uint16_t i = 0;
    for (auto kmer_pair : kmer_pairs)
    {
        if (filter_cross_dimerization(kmer_pair.get_kmer_ID1(), kmer_pair.get_kmer_ID2()))
            mask.set(i);
        ++i;
    }
}

// post-filter candidates fulfilling chemical constraints by their relative frequency
/*void post_frequency_filter(TKmerLocations kmer_locations, TSeqNo occurrence_freq)
{

}*/

// filter k-mers by frequency and chemical properties
void pre_filter_main(io_cfg_type const & io_cfg, primer_cfg_type const & primer_cfg, TKLocations & locations, TKmerLocations & kmer_locations)
{
    using TSeqNo = typename seqan::Value<typename TLocations::key_type, 1>::Type;

    // scale to be lower frequency bound for filters
    TSeqNo cutoff = primer_cfg.cutoff;
    // continue here
    std::cout << "INFO: Cut-off frequency = " << cutoff << std::endl;
    // frequency filter and sequence fetching
    frequency_filter(io_cfg, locations, kmer_locations, cutoff);
    std::cout << "INFO: kmers after frequency cutoff = " << kmer_locations.size() << std::endl;
    chemical_filter_single(primer_cfg, kmer_locations);
    std::cout << "INFO: kmers after chemical filtering = " << kmer_locations.size() << std::endl;
}

// combine helper, forward window until both iterators point to same reference ID or end
void fast_forward(TKmerLocation::TLocationVec const & locations1, TKmerLocations::value_type::const_iterator & it1_loc_start, TKmerLocation::TLocationVec const & locations2, TKmerLocations::value_type::const_iterator & it2_loc_start)
{
    //std::cout << "enter fast_forward ..." << std::endl;
    while (it1_loc_start != locations1.end() && it2_loc_start != locations2.end() &&
            seqan::getValueI1<TSeqNo, TSeqPos>(*it1_loc_start) != seqan::getValueI1<TSeqNo, TSeqPos>(*it2_loc_start))
    {
        // speedup: use lower_bound/upper_bound
        while (it1_loc_start != locations1.end() && seqan::getValueI1<TSeqNo, TSeqPos>(*it1_loc_start) < seqan::getValueI1<TSeqNo, TSeqPos>(*it2_loc_start))
            ++it1_loc_start;
        while (it2_loc_start != locations2.end() && seqan::getValueI1<TSeqNo, TSeqPos>(*it2_loc_start) < seqan::getValueI1<TSeqNo, TSeqPos>(*it1_loc_start))
            ++it2_loc_start;
    }
    //std::cout << "... leaving ff" << std::endl;
}

/*
 * Helper for combine.
 *
 */

/* Combine based on suitable location distances s.t. transcript length is in permitted range.
 * Chemical suitability will be tested by a different function. First position indicates,
 * that the k-mer corresponds to a forward primer, and second position indicates reverse
 * primer, i.e. (k1, k2) != (k2, k1).
 */
void combine(primer_cfg_type const & primer_cfg, TKmerLocations const & kmer_locations, TKmerPairs & kmer_pairs)
{
    using it_loc_type = TKmerLocations::value_type::const_iterator;
    it_loc_type it1_loc_start, it1_loc_aux;
    it_loc_type it2_loc_start, it2_loc_aux;
    TKmer kmer1, kmer2;
    TKmerID kmer_ID1, kmer_ID2;
    TKmerLength K1, K2;
    for (auto it1 = kmer_locations.begin(); it1 != kmer_locations.end() && it1 != kmer_locations.end()-1; ++it1)
    {
        K1 = (*it1).get_K();
        for (auto it2 = it1+1; it2 != kmer_locations.end(); ++it2)
        {
            K2 = (*it2).get_K();
            kmer_ID1 = (*it1).get_kmer_ID();
            kmer_ID2 = (*it2).get_kmer_ID();
            assert(kmer_ID1 && kmer_ID2);

            // continue with next combination if kmer sequences do not pass cross-dimerization filter
            if (!filter_cross_dimerization(kmer_ID1, kmer_ID2))
                continue;
            // iterator to start position of current location for k-mer 1
            it1_loc_start = (*it1).locations.begin();
            // iterator to start position of current location for k-mer 2
            it2_loc_start = (*it2).locations.begin();
            // forward iterators to correspond to refer to same sequence ID or end
            fast_forward((*it1).locations, it1_loc_start, (*it2).locations, it2_loc_start);
            // no common reference sequences => forward kmer iterators
            if (it1_loc_start == (*it1).locations.end() || it2_loc_start == (*it2).locations.end())
                continue;

            // foward iterator for k-mer 1 on same sequence (loc)
            it1_loc_aux = it1_loc_start;
            // foward iterator for k-mer 2 on same sequence (loc)
            it2_loc_aux = it2_loc_start;
            // invariant after entering this loop: loc(it1) == loc(it2)
            auto seq_ID = seqan::getValueI1<TSeqNo, TSeqPos>(*it1_loc_start);
            while (it1_loc_start != (*it1).locations.end() && it1_loc_aux != (*it1).locations.end() && it2_loc_start != (*it2).locations.end() && seq_ID == seqan::getValueI1<TSeqNo, TSeqPos>(*it2_loc_aux))
            {
                // valid combination?
                auto pos_kmer1 = seqan::getValueI2<TSeqNo, TSeqPos>(*it1_loc_aux);
                if (it2_loc_aux != (*it2).locations.end())
                {
                    auto pos_kmer2 = seqan::getValueI2<TSeqNo, TSeqPos>(*it2_loc_aux);
                    auto pos_delta = (pos_kmer1 < pos_kmer2) ? pos_kmer2 - pos_kmer1 - K1: pos_kmer1 - pos_kmer2 - K2;
                    if (pos_delta >= primer_cfg.get_transcript_range().first && pos_delta <= primer_cfg.get_transcript_range().second)
                    {
                        TKmerID kmer_fwd_new = (pos_kmer1 < pos_kmer2) ? (*it1).get_kmer_ID() : (*it2).get_kmer_ID();
                        TKmerID kmer_rev_new = (pos_kmer1 < pos_kmer2) ? (*it2).get_kmer_ID() : (*it1).get_kmer_ID();
                        auto pair_location = std::make_tuple(seq_ID, std::min<TSeqPos>(pos_kmer1, pos_kmer2), std::max<TSeqPos>(pos_kmer1, pos_kmer2));

                        // extend location vector if pair combinations already in result
                        if (kmer_pairs.size() && kmer_pairs.back().get_kmer_ID1() == kmer_fwd_new && kmer_pairs.back().get_kmer_ID2() == kmer_rev_new)
                        {
                            kmer_pairs[kmer_pairs.size()-1].pair_locations.push_back(pair_location);
                        }
                        else
                        {
                            TKmerPair pair{kmer_fwd_new, kmer_rev_new, abs(primer_melt_wallace(kmer_fwd_new) - primer_melt_wallace(kmer_rev_new)), pair_location};
                            kmer_pairs.push_back(pair);
                        }
                    }
                    ++it2_loc_aux;
                }
                // all combinations tested for second k-mer
                if (it2_loc_aux == (*it2).locations.end())
                {
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
            }
        }
    }
}

/*void post_filter_main(primer_cfg_type const & primer_cfg, TKmerLocations & kmer_locations, TKmerPairs & pairs)
{

}*/

}  // namespace priset
