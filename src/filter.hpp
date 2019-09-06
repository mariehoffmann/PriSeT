#pragma once

#include <algorithm>
#include <bitset>
#include <cmath>
#include <fstream>
#include <functional>
#include <string>
#include <unordered_set>

#include <sdsl/bit_vectors.hpp>

#include "../submodules/genmap/src/common.hpp"
#include "../submodules/genmap/src/genmap_helper.hpp"

//#include "../submodules/sdsl/include/sdsl/bit_vectors.hpp"

#include "primer_cfg_type.hpp"
#include "types.hpp"
#include "utilities.hpp"

namespace priset
{
/*
 * Filter kmers based on their chemical properties regardless of their pairing.
 * Constraints that are checked: melting tempaerature, CG content
 */
bool chemical_filter_single(primer_cfg_type const & primer_cfg, TKmerID const & kmerID)
{
    assert(kmerID > 0);

    float Tm = get_Tm(primer_cfg, kmerID);

    // Filter by melting temperature
    if (Tm >= primer_cfg.get_min_Tm() && Tm <= primer_cfg.get_max_Tm())
    {
        // Filter for CG content.
        if (filter_CG(primer_cfg, kmerID))
        {
            // Filter if Gibb's free energy is below -6 kcal/mol
            if (filter_self_dimerization(kmerID))
            {
                if (filter_repeats_runs(kmerID))
                    return true;
            }

        }
    }

/*
    for (TKmerLocations::size_type i = 0; i < kmer_locations.size(); ++i)
    {

        // TODO: optimize - drop kmer when computing Tm in sequence lookup fct
        kmerID = kmer_locations[i].get_kmer_ID1();
        //std::cout << "kmerID = " << kmerID << std::endl;
        auto Tm = primer_melt_wallace(kmerID);
    //    std::cout << "Tm = " << Tm << std::endl;
        // filter by melting temperature
        if (Tm >= Tm_min && Tm <= Tm_max)
        {
            // filter by CG content
    //        std::cout << "call filter_CG ...\n";
            if (filter_CG(primer_cfg, kmerID))
            {

                // Filter if Gibb's free energy is below -6 kcal/mol
                if (filter_self_dimerization(kmerID))
                {
                    if (filter_repeats_runs(kmerID))
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
*/
    return false;
}


// TODO: kmer_location is the transformed to reference bit vector
// pre-filter and sequence fetch
// 1. filter candidates by number of occurences only independent of their chemical suitability
// 2. fetch sequence and check chemical constraints that need to hold for a single primer
/*void frequency_filter(priset::io_cfg_type const & io_cfg, TKLocations & locations, TKmerLocations & kmer_locations, TSeqNo const cutoff)
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
*/
// helper for dna sequence extraction
/*template<typename TText>
TKmerID encode_wrapper(TText const & text, TSeqNo const seqNo, TSeqPos const seqPos, TKmerLength const K)
{
    seqan::DnaString seq = seqan::valueById(text, seqNo);
    auto const & kmer_str = seqan::infixWithLength(seq, seqPos, K);
    return dna_encoder(kmer_str);
};*/

// Filter of single kmers and transform of references to bit vectors.
void filter_and_transform(io_cfg_type const & io_cfg, primer_cfg_type const & primer_cfg, TKLocations & locations, TReferences & references, TKmerIDs & kmerIDs, TSeqNoMap & seqNoMap, TSeqNo const cutoff, TKmerCounts & stats)
{
    //filter_stats stats{};
    std::cout << "Enter filter_and_transform ...\n";
    // uniqueness indirectly preserved by (SeqNo, SeqPos) if list sorted lexicographically
    assert(length(locations));

    references.clear();
    kmerIDs.clear();

    // load corpus for dna to 64 bit conversion
    seqan::StringSet<seqan::DnaString, seqan::Owner<seqan::ConcatDirect<>>> text;
    fs::path text_path = io_cfg.get_index_txt_path();
    std::cout << "text_path = " << text_path << std::endl;
    seqan::open(text, text_path.string().c_str(), seqan::OPEN_RDONLY);
    std::cout << "text length = " << seqan::length(text) << std::endl;
    for (auto ss : text)
        std::cout << ss << ", ";
    std::cout << std::endl;

    // (i) collect distinct sequence identifiers and maximal position of kmer occurences
    // to have a compressed representation.
    std::cout << "Step 1\n";
    std::vector<TSeqPos> seqNo2maxPos(1 << 10);
    for (typename TKLocations::const_iterator it = locations.begin(); it != locations.end(); ++it)
    {
    //    std::cout << "iterate through locations\n";
        // cleanup in mapper may lead to undercounting kmer occurrences
        if ((it->second).first.size() < cutoff)
            continue;

        const auto & [seqNo, seqPos, K] = (it->first);
        // use symmetry and lexicographical ordering of locations to skip already seen ones
        // TODO: is this already shortcut fm mapper?
        if (it->second.first.size() && (seqan::getValueI1<TSeqNo, TSeqPos>(it->second.first[0]) < seqNo || (seqan::getValueI1<TSeqNo, TSeqPos>(it->second.first[0]) == seqNo &&  seqan::getValueI2<TSeqNo, TSeqPos>(it->second.first[0]) < seqPos)))
            continue;
        // update largest kmer position
        if (seqNo >= seqNo2maxPos.size())
            seqNo2maxPos.resize(seqNo + 1);
        seqNo2maxPos[seqNo] = std::max(seqNo2maxPos.at(seqNo), seqPos);
        // and also for other occurences
        for (std::vector<TLocation>::const_iterator it_loc_fwd = it->second.first.begin(); it_loc_fwd != it->second.first.end(); ++it_loc_fwd)
        {
            TSeqNo seqNo_fwd = seqan::getValueI1<TSeqNo, TSeqPos>(*it_loc_fwd);
            TSeqPos seqPos_fwd = seqan::getValueI2<TSeqNo, TSeqPos>(*it_loc_fwd);
            if (seqNo_fwd >= seqNo2maxPos.size())
                seqNo2maxPos.resize(seqNo_fwd + 1);
            seqNo2maxPos[seqNo_fwd] = std::max(seqNo2maxPos.at(seqNo_fwd), seqPos_fwd);
        }
    }

    // (ii) Create bit vectors in the length of largest kmer occurences, and set
    // bits for kmer occurrences. Collect also kmer lengths encoded in final kmer
    // code of type uint64_t. Final sequence lookup and encoding is done in next step.
    std::cout << "max positions for sequences: \nseqNo\tmax(seqPos)\n";
    for (unsigned i = 0; i < seqNo2maxPos.size(); ++i)
        if (seqNo2maxPos[i])
            std::cout << i << "\t->\t" << seqNo2maxPos[i] << std::endl;
    seqNoMap.clear();

    unsigned i_cx = 0;
    // reserve space for bit vectors
    for (unsigned i = 0; i < seqNo2maxPos.size(); ++i)
    {
        if (!seqNo2maxPos[i])
            continue;
        std::cout << "seqNoMap, insert: " << i << " -> " << seqNo2maxPos[i] << std::endl;
        seqNoMap[i] = i_cx++;
        sdsl::bit_vector bv(seqNo2maxPos[i] + 1, 0);
        references.push_back(bv);
    }

    // set bits for kmers in references and collect Ks
    // bit 0: set if kmer of length primer_length_min found ...
    // bit primer_length_max-primer_length_min: set if kmer of length primer_length_max found
    // => maximal encodable length difference is 16!
    std::unordered_map<uint64_t, uint64_t> loc_and_ks;
    kmerIDs.clear();
    kmerIDs.resize(references.size());
    std::vector<uint32_t> debug_drop_kmer_repeats(references.size(), 0);
    //std::vector<std::vector<std::pair<TSeqPos, TKmerID>>> pos_kmers(references.size());
    for (typename TKLocations::const_iterator it = locations.begin(); it != locations.end(); ++it)
    {
        // cleanup in mapper may lead to undercounting kmer occurrences
        if ((it->second).first.size() < cutoff)
            continue;
        const auto & [seqNo, seqPos, K] = (it->first);
        std::cout << "CURRENT K = " << K << std::endl;
        // use symmetry and lexicographical ordering of locations to skip already seen ones
        // TODO: is this already shortcut fm mapper?
        if (it->second.first.size() && (seqan::getValueI1<TSeqNo, TSeqPos>(it->second.first[0]) < seqNo ||
            (seqan::getValueI1<TSeqNo, TSeqPos>(it->second.first[0]) == seqNo &&
            seqan::getValueI2<TSeqNo, TSeqPos>(it->second.first[0]) < seqPos)))
            continue;

        // insert now unique occurrences listed in first vector (forward mappings)
        TSeqNo seqNo_prev{0};
        TSeqPos seqPos_prev{0};
        for (std::vector<TLocation>::const_iterator it_loc_fwd = it->second.first.begin(); it_loc_fwd != it->second.first.end(); ++it_loc_fwd)
        {
            TSeqNo seqNo = seqan::getValueI1<TSeqNo, TSeqPos>(*it_loc_fwd);
            TSeqPos seqPos = seqan::getValueI2<TSeqNo, TSeqPos>(*it_loc_fwd);
            auto seqNo_cx = seqNoMap.at(seqNo); // compressed sequence id
            // continue and delete previous bit if same kmer occurs within 400 bp
            if (it_loc_fwd > it->second.first.begin() && seqNo_prev == seqNo && seqPos_prev + primer_cfg.same_kmer_distance < seqPos)
            {
                // undo previous kmer bit in reference
                references[seqNo_cx][seqPos_prev] = 0;
                seqPos_prev = seqPos;
                debug_drop_kmer_repeats[seqNo_cx]++;
                continue;
            }

            references[seqNo_cx][seqPos] = 1;
            auto K_store_offset = K - primer_cfg.primer_min_length;
            if (K_store_offset > 12)
                throw std::invalid_argument("ERROR: kmer length difference exceeds 12 + primer_min_length - 1 bp!");
            uint64_t loc_key = location_encode(seqNo, seqPos);
            uint64_t loc_val = 1ULL << (63ULL - K_store_offset);
            // locations filled for K_min to K_max
            if (loc_and_ks.find(loc_key) == loc_and_ks.end())
                loc_and_ks.insert({loc_key, loc_val});
            else // update
                loc_and_ks[loc_key] += loc_val;
            seqNo_prev = seqNo;
            seqPos_prev = seqPos;
        }
    }

    std::cout << "INFO: dropped kmers due to close repetetion: ";
    for (auto ctr : debug_drop_kmer_repeats)
        std::cout << ctr << " " << std::endl;

    // (iii) lookup kmer sequences, filter and encode as 64 bit integers.
    // inverse reference ID map
    // todo: use one bit vector instead with rank support
    std::unordered_map<TSeqNo, TSeqNo> seqNoMap_inv;
    for (auto it = seqNoMap.begin(); it != seqNoMap.end(); ++it)
        seqNoMap_inv.insert({it->second, it->first});
    kmerIDs.resize(references.size());
    uint64_t const kmer_length_mask = ~((1ULL << 52ULL) - 1ULL);
    for (unsigned i = 0; i < references.size(); ++i)
    {
        sdsl::rank_support_v5<1> r1s(&references[i]); // check if once initialized modification of references[i] does not invalidate rank support
        sdsl::select_support_mcl<1,1> s1s(&references[i]);
        for (unsigned r = r1s.rank(references[i].size()); r > 0; --r)
        {
            TSeqNo seqNo = seqNoMap_inv[i];
            TSeqPos seqPos = s1s.select(r);
            uint64_t loc_key = location_encode(seqNo, seqPos);
            std::cout << "seqNo and seqpos as key = " << loc_key << std::endl;
            if (loc_and_ks.find(loc_key) == loc_and_ks.end())
                throw std::invalid_argument("ERROR: " + std::to_string(loc_key) + " not in loc_and_ks dictionary.");

            TKmerID kmerID_head = loc_and_ks[loc_key]; // tail not filled yet & kmer_length_mask;
            TKmerID k_pattern = kmerID_head >> 52ULL;
            TKmerID k_pattern_cpy{k_pattern};
            std::cout << "k_pattern = " << k_pattern << std::endl;
            //while(!k_pattern_cpy)
            //    std::cout << ((k_pattern  & 1) ? "1" : "0");
            // TODO: identify highest set bit in head
            TKmerLength k = primer_cfg.primer_min_length + primer_cfg.pattern_word_length - 1;
            while (!(k_pattern_cpy & 1) && k--)
                k_pattern_cpy >>= 1;

            std::cout << "MSG: maximal identified k = " << k << std::endl;
            // lookup sequence in corpus and encode
            seqan::DnaString seq = seqan::valueById(text, seqNo);
            TSeq const & kmer_str = seqan::infixWithLength(seq, seqPos, k);
            TKmerID kmerID = dna_encoder(kmer_str);
            TKmerID kmerID_trim = kmerID;
            kmerID |= kmerID_head;
            TKmerID trim_offset = 0;
            k_pattern_cpy = k_pattern; // reset copy
            while (!k_pattern)
            {
                if (k_pattern & 1)
                {
                    // erase bit in head
                    if (!chemical_filter_single(primer_cfg, kmerID_trim))
                    {
                        kmerID ^= 1 << (52ULL + trim_offset);
                    }
                    else
                        stats[KMER_COUNTS::FILTER1_CNT]++;
                }
                kmerID_trim >>= 2;
                k_pattern >>= 1;
                ++trim_offset;
            }
            // do not store Kmer and reset bit if for no length the filter was passed
            if (!(kmer_length_mask & kmerID))
                references[i][seqPos] = 0;
            else
                kmerIDs[i].push_front(kmerID);
        }
    }
    // TODO: delete references with no more bits, or too inefficient w.r.t. possible space gain?

    for (auto kmerID_list : kmerIDs)
    {
        for (auto kmerID : kmerID_list)
            if (!kmerID)
            {
                std::cout << "ERROR: kmerID = 0" << std::endl;
                exit(0);
            }
    }
    locations.clear();
}

// check cross-dimerization.
// Note: No length bit shall be set in head, code corresponds to true length.
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
void pre_filter_main(io_cfg_type const & io_cfg, primer_cfg_type const & primer_cfg, TKLocations & locations, TReferences & references, TKmerIDs & kmerIDs, TSeqNoMap & seqNoMap, TKmerCounts & kmerCounts)
{
    using TSeqNo = typename seqan::Value<typename TLocations::key_type, 1>::Type;

    // scale to be lower frequency bound for filters
    TSeqNo cutoff = primer_cfg.cutoff;
    // continue here
    std::cout << "INFO: Cut-off frequency = " << cutoff << std::endl;
    // frequency filter and sequence fetching
    filter_and_transform(io_cfg, primer_cfg, locations, references, kmerIDs, seqNoMap, primer_cfg.cutoff, kmerCounts);
    std::unordered_set<TKmerID> kmer_set;
    unique_kmers(kmerIDs, kmer_set);
    std::cout << "INFO: kmers after frequency cutoff and chemical filtering = " << kmer_set.size() << std::endl;
}

/*
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
*/


/* Combine based on suitable location distances s.t. transcript length is in permitted range.
 * Chemical suitability will be tested by a different function. First position indicates,
 * that the k-mer corresponds to a forward primer, and second position indicates reverse
 * primer, i.e. (k1, k2) != (k2, k1).
 */
 /*
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
*/

// primer_cfg_type const & primer_cfg, TKmerLocations const & kmer_locations, TKmerPairs & kmer_pairs
void combine2(primer_cfg_type const & primer_cfg, TReferences const & references, TKmerIDs const & kmerIDs, TPairs & pairs, TKmerCounts & stats)
{
    pairs.clear();
    pairs.resize(kmerIDs.size());
    uint64_t offset_max = primer_cfg.get_transcript_range().second;

    for (uint64_t i = 0; i < references.size(); ++i)
    {
        std::cout << "reference ID = " << i << std::endl;
        //TReference reference = references.at(i);
        //auto kmerID_list = kmerIDs.at(i);

        sdsl::bit_vector reference;
        sdsl::util::assign(reference, references.at(i));
        sdsl::rank_support_v5<1, 1> r1s(&reference); // replace after bugfix with
        sdsl::select_support_mcl<1> s1s(&reference);

        std::cout << "iterate over " << r1s.rank(reference.size()) << " kmer IDs\n";
        for (uint64_t r_fwd = 1; r_fwd <= r1s.rank(reference.size()); ++r_fwd)
        {
            // text position of r-th k-mer
            uint64_t idx_fwd = s1s.select(r_fwd);
            std::cout << "text index of r_fwd = " << r_fwd << " is " << idx_fwd << std::endl;
            TKmerID kmerID_fwd = kmerIDs[i][r_fwd-1];
            if (!(kmerID_fwd >> primer_cfg.kmer_word_length))
            {
                std::cout << "ERROR: k length pattern is zero\n";
                exit(-1);
            }
            std::cout << "INFO: kurrent kmerID_fwd = " << kmerID_fwd << std::endl;

            TKmerID mask_fwd = 1ULL << 63ULL;
            while (mask_fwd >= (1ULL << primer_cfg.kmer_word_length))
            {
                std::cout << "\tINFO: current length bit mask = " << print_bits(mask_fwd);
                if (mask_fwd & kmerID_fwd)
                {
                    // window start position (inclusive)
                    auto k_fwd = 63 - log2_asm(mask_fwd) + primer_cfg.primer_min_length;
                    uint64_t w_begin = idx_fwd + k_fwd + primer_cfg.get_transcript_range().first + 1;
                    std::cout << "\tINFO: current length pattern in kmerID_fwd corresponding to k_fwd = " << k_fwd << std::endl;
                    // Note: references are truncated to the last 1-bit, hence last window start position is reference.size()-1
                    if (w_begin > reference.size() - 1)
                        break;
                    // window end position (exclusive)
                    uint64_t w_end = std::min(reference.size(), idx_fwd + k_fwd + offset_max + 1);
                    std::cout << "\tINFO: search window = [" << w_begin << ", " << w_end << "], rank(w_begin) = " << r1s.rank(w_begin) << ", rank(w_end) = " << r1s.rank(w_end) << " corresponding to " << (r1s.rank(w_end) - r1s.rank(w_begin)) << " candidate pairing mates\n";
                    // iterate through kmers in reference sequence window [w_begin : w_end]
                    for (uint64_t r_rev = r1s.rank(w_begin); r_rev <= r1s.rank(w_end); ++r_rev)
                    {
                        // iterate through encoded lengths
                        TKmerID kmerID_rev = kmerIDs.at(i).at(r_rev - 1);
                        TKmerID mask_rev = 1ULL << 63ULL;
                        TCombinePattern cp{};
                        std::cout << "\t\tINFO: current kmerID_rev = " << kmerID_rev << std::endl;
                        while (mask_rev >= (1ULL << primer_cfg.kmer_word_length))
                        {
                            std::cout << "\t\tINFO: current length bit mask = ";
                            print_bits(mask_fwd);
                            if (mask_rev & kmerID_rev)
                            {
                                std::cout << "\t\tINFO: current length pattern in kmerID_fwd corresponding to k_fwd = " << log2_asm(mask_rev) << std::endl;
                                // TODO: add more filter here
                                std::cout << "\t\tINFO: compute Tm_delta(" << kmerID_fwd << ", " << (63 - log2_asm(mask_fwd)+primer_cfg.primer_min_length) << ", " << kmerID_rev << ", " << (63 - log2_asm(mask_rev) + primer_cfg.primer_min_length) << ") = " << Tm_delta(kmerID_fwd, mask_fwd, kmerID_rev, mask_rev) << std::endl;
                                if (Tm_delta(kmerID_fwd, mask_fwd, kmerID_rev, mask_rev) <= primer_cfg.primer_melt_diff)
                                {
                                    std::cout << "\t\tINFO: Tm_delta in range, set combination bit\n";
                                    cp.set(mask_fwd, mask_rev, primer_cfg.pattern_word_length);
                                    ++stats[KMER_COUNTS::COMBINER_CNT];
                                }
                                else
                                    std::cout << "\t\tINFO: Tm_delta outside of range, no bit setting\n";
                            }
                            mask_rev >>= 1ULL;
                        }
                        if (cp.data[0] | cp.data[1])
                            pairs[i].push_back(TPair{r_fwd, r_rev, cp});
                    }
                }
                mask_fwd >>= 1ULL;
            }
        }
        std::cout << "next reference\n";
    }
    std::cout << "leaving combine2\n";
}

/*
void post_filter_main(primer_cfg_type const & primer_cfg, TKmerIDs & kmerIDs, TPairs & pairs)
{

}*/

}  // namespace priset
