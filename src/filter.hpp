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

#include "combine_types.hpp"
#include "primer_cfg_type.hpp"
#include "types.hpp"
#include "utilities.hpp"

namespace priset
{
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
    //std::cout << "Enter filter_and_transform ...\n";
    // uniqueness indirectly preserved by (SeqNo, SeqPos) if list sorted lexicographically
    assert(length(locations));

    references.clear();
    kmerIDs.clear();

    // load corpus for dna to 64 bit conversion
    seqan::StringSet<seqan::DnaString, seqan::Owner<seqan::ConcatDirect<>>> text;
    fs::path text_path = io_cfg.get_index_txt_path();
    //std::cout << "text_path = " << text_path << std::endl;
    seqan::open(text, text_path.string().c_str(), seqan::OPEN_RDONLY);
    std::cout << "text length = " << seqan::length(text) << std::endl;
    for (auto ss : text)
        std::cout << ss << ", ";
    std::cout << std::endl;

    // (i) collect distinct sequence identifiers and maximal position of kmer occurences
    // to have a compressed representation.
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
            if (it_loc_fwd > it->second.first.begin() && seqNo_prev == seqNo && seqPos_prev + TRAP_DIST < seqPos)
            {
                // undo previous kmer bit in reference
                references[seqNo_cx][seqPos_prev] = 0;
                seqPos_prev = seqPos;
                debug_drop_kmer_repeats[seqNo_cx]++;
                continue;
            }

            references[seqNo_cx][seqPos] = 1;
            std::cout << "DEBUG: K = " << K << ", PRIMER_MIN_LEN = " << PRIMER_MIN_LEN << ", PRIMER_MAX_LEN = " << PRIMER_MAX_LEN << std::endl;
            if (uint64_t(K) > PRIMER_MAX_LEN)
                throw std::invalid_argument("ERROR: kmer length difference exceeds 12 + primer_min_length - 1 bp!");
            uint64_t loc_key = location_encode(seqNo, seqPos);
            uint64_t loc_val = ONE_LSHIFT_63 >> (K - PRIMER_MIN_LEN);
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
    std::cout << "here1\n";
    std::unordered_map<TSeqNo, TSeqNo> seqNoMap_inv;
    for (auto it = seqNoMap.begin(); it != seqNoMap.end(); ++it)
        seqNoMap_inv.insert({it->second, it->first});
    kmerIDs.resize(references.size());
    std::cout << "here2\n";
    uint64_t const kmer_length_mask = ~((1ULL << KMER_SIZE) - 1ULL);
    for (unsigned i = 0; i < references.size(); ++i)
    {
        std::cout << "here3 loop start\n";
        sdsl::rank_support_v5<1> r1s(&references[i]); // check if once initialized modification of references[i] does not invalidate rank support
        sdsl::select_support_mcl<1,1> s1s(&references[i]);
        for (unsigned r = r1s.rank(references[i].size()); r > 0; --r)
        {
            TSeqNo seqNo = seqNoMap_inv[i];
            TSeqPos seqPos = s1s.select(r);
            uint64_t loc_key = location_encode(seqNo, seqPos);
            if (loc_and_ks.find(loc_key) == loc_and_ks.end())
                throw std::invalid_argument("ERROR: " + std::to_string(loc_key) + " not in loc_and_ks dictionary.");

            TKmerID kmerID_prefix = loc_and_ks[loc_key]; // tail not filled yet & KMER_SIZEgth_mask;
            //TKmerID k_pattern = kmerID_prefix >> (WORD_SIZE - LEN_MASK_SIZE);
            //TKmerID k_pattern_cpy{k_pattern};
            // identify lowest set bit in head
            TKmerLength k_max = PRIMER_MAX_LEN - ffsll(kmerID_prefix >> 54) + 1;

            std::cout << "here4, seqNo = " << seqNo << std::endl;
            // lookup sequence in corpus and encode
            seqan::DnaString seq = seqan::valueById(text, seqNo);
            std::cout << "here4a: seq = " << seq << ", seqPos = " << seqPos << ", k_max = " << k_max << std::endl;

            TSeq const & kmer_str = seqan::infixWithLength(seq, seqPos, k_max);
            std::cout << "here4b, kmer_str = " << kmer_str << std::endl;
            TKmerID kmerID = kmerID_prefix + dna_encoder(kmer_str);

            std::cout << "here5\n";
            // erase those length bit in prefix corresponding to kmers not passing the filter
            chemical_filter_single_pass(primer_cfg, kmerID);
            std::cout << "here6\n";
            // do not store Kmer and reset bit in reference
            if (!(MASK_SELECTOR & kmerID))
                references[i][seqPos] = 0;
            else
            {
                stats[KMER_COUNTS::FILTER1_CNT] += __builtin_popcountll(kmerID >> 54);
                kmerIDs[i].push_front(kmerID);
            }
        }
    }
    // TODO: delete references with no more bits, or too inefficient w.r.t. possible space gain?
    std::cout << "here7\n";
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
    std::cout << "leaving filter and transform\n";
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
}

/* Combine based on suitable location distances s.t. transcript length is in permitted range.
 * Chemical suitability will be tested by a different function. First position indicates,
 * that the k-mer corresponds to a forward primer, and second position indicates reverse
 * primer, i.e. (k1, k2) != (k2, k1).
 */
// primer_cfg_type const & primer_cfg, TKmerLocations const & kmer_locations, TKmerPairs & kmer_pairs
// TODO: add concept requring outer and inner container to provide begin/end/push_back
template<typename TPairList>
void combine(primer_cfg_type const & primer_cfg, TReferences const & references, TKmerIDs const & kmerIDs, TPairList & pairs, TKmerCounts & stats)
{
    pairs.clear();
    for (uint64_t i = 0; i < references.size(); ++i)
    {
        std::cout << "DEBUG: reference ID = " << i << std::endl;

        sdsl::bit_vector reference;
        sdsl::util::assign(reference, references.at(i));
        sdsl::rank_support_v5<1, 1> r1s(&references.at(i)); // replace after bugfix with
        sdsl::select_support_mcl<1> s1s(&reference);

        std::cout << "DEBUG: iterate over " << r1s.rank(reference.size()) << " kmer IDs\n";
        for (uint64_t r_fwd = 1; r_fwd < r1s.rank(reference.size()); ++r_fwd)
        {
            // text position of r-th k-mer
            uint64_t idx_fwd = s1s.select(r_fwd);
            std::cout << "DEBUG: text index of r_fwd = " << r_fwd << " is " << idx_fwd << std::endl;
            TKmerID const kmerID_fwd = kmerIDs[i][r_fwd - 1];
            if (!(kmerID_fwd >> KMER_SIZE))
            {
                std::cout << "ERROR: k length pattern is zero\n";
                exit(-1);
            }
            // minimal window start position for pairing kmer
            uint64_t w_begin = idx_fwd + PRIMER_MIN_LEN + TRANSCRIPT_MIN_LEN;

            // maximal window end position (exclusive) for pairing kmer
            uint64_t w_end = std::min(reference.size(), idx_fwd + PRIMER_MAX_LEN + TRANSCRIPT_MAX_LEN + 1);

            // iterate through kmers in reference sequence window [w_begin : w_end]
            // note that w_begin/end are updated due to varying kmer length of same kmerID
            for (uint64_t r_rev = r1s.rank(w_begin) + 1; r_rev <= r1s.rank(w_end); ++r_rev)
            {
                // new kmer combination, store possible combinations
                TCombinePattern<TKmerID, TKmerLength> cp; // data set to 0?
                // reset mask to highest bit (= PRIMER_MIN_LEN position)
                uint64_t mask_fwd = ONE_LSHIFT_63;
                while ((((mask_fwd << 1) - 1) & kmerID_fwd) >> (WORD_SIZE - LEN_MASK_SIZE))
                {
                    if (mask_fwd & kmerID_fwd)
                    {
                        std::cout << "INFO: current k of kmerID_fwd = " << __builtin_clzl(mask_fwd)+PRIMER_MIN_LEN << std::endl;
                        // window start position (inclusive)

                        std::cout << "\tDEBUG: search window = [" << w_begin << ", " << w_end << "], rank(w_begin) = " << r1s.rank(w_begin) << ", rank(w_end) = " << r1s.rank(w_end) << " corresponding to " << (r1s.rank(w_end) - r1s.rank(w_begin)) << " candidate pairing mates\n";

                        // iterate through encoded lengths
                        TKmerID const kmerID_rev = kmerIDs.at(i).at(r_rev - 1);
                        uint64_t mask_rev = ONE_LSHIFT_63;
                        // as long there are bits in length mask part of kmer_rev we continue shifting mask_rev
                        while ((((mask_rev << 1) - 1) & kmerID_rev) >> (WORD_SIZE - LEN_MASK_SIZE))
                        {
                            if (mask_rev & kmerID_rev)
                            {
                                std::cout << "\t\tDEBUG: current length pattern in kmerID_rev corresponding to k_rev = " << (__builtin_clzl(mask_rev) + PRIMER_MIN_LEN) << std::endl;
                                // TODO: add more filter here
                                std::cout << "\t\tINFO: compute dTm(" << dna_decoder(kmerID_fwd, mask_fwd) << ", " << dna_decoder(kmerID_rev, mask_rev) << ") = ";
                                std::cout << dTm(kmerID_fwd, mask_fwd, kmerID_rev, mask_rev) << std::endl;
                                if (dTm(kmerID_fwd, mask_fwd, kmerID_rev, mask_rev) <= primer_cfg.primer_melt_diff)
                                {
                                    std::cout << "\t\tINFO: dTm in range, set combination bit\n";
                                    cp.set(mask_fwd, mask_rev);
                                    ++stats[KMER_COUNTS::COMBINER_CNT];
                                }
                                else
                                    std::cout << "\t\tINFO: dTm outside of range, no bit setting\n";
                            }
                            mask_rev >>= 1ULL; // does not affect search window, since starting position is fixed
                        } // length mask_rev
                    }
                    // update mask and dependent search window boundaries
                    mask_fwd >>= 1ULL;
                    ++w_begin;
                    // Note: references are truncated to the last 1-bit, hence last window start position is reference.size()-1
                    // window end position (exclusive)
                    // continue with next kmer combination since we are outside of search window
                    if ((w_begin > reference.size() - 1) || (w_end < w_begin + TRANSCRIPT_MIN_LEN))
                        break;
                } // length mask_fwd
                if (cp.is_set())
                {
                    std::cout << "\t\tINFO: push back combination pattern\n";
                    pairs.push_back(TPair<TCombinePattern<TKmerID, TKmerLength>>{i, r_fwd, r_rev, cp});
                }
                else
                    std::cout << "\t\tINFO: did NOT push back combination pattern\n";
            } // kmerID rev
        } // kmerID fwd
        std::cout << "next reference\n";
    }
    std::cout << "leaving combine\n";
}

}  // namespace priset
