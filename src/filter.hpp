#pragma once

#include <algorithm>
#include <bitset>
#include <cmath>
#include <fstream>
#include <functional>
#include <string>
#include <unordered_set>

#include "../submodules/genmap/src/common.hpp"
#include "../submodules/genmap/src/genmap_helper.hpp"
#include "../submodules/sdsl-lite/include/sdsl/bit_vectors.hpp"

#include "primer_cfg_type.hpp"
#include "types.hpp"
#include "utilities.hpp"

namespace priset
{

// Filter of single kmers and transform of references to bit vectors.
void transform_and_filter(io_cfg_type const & io_cfg, TKLocations & locations, TReferences & references, TSeqNoMap & seqNoMap, TKmerIDs & kmerIDs, TKmerCounts * kmerCounts = nullptr)
{
    // uniqueness indirectly preserved by (SeqNo, SeqPos) if list sorted lexicographically
    assert(length(locations));

    references.clear();
    kmerIDs.clear();

    // load corpus for dna to 64 bit conversion
    seqan::StringSet<seqan::DnaString, seqan::Owner<seqan::ConcatDirect<>>> text;
    fs::path text_path = io_cfg.get_index_txt_path();
    seqan::open(text, text_path.string().c_str(), seqan::OPEN_RDONLY);
    if (!seqan::lengthSum(text))
        throw std::length_error("Reference text size is 0!");
    // (i) collect distinct sequence identifiers and maximal position of kmer occurences
    // to have a compressed representation.
    // std::map<TSeqNo, TSeqPos> seqNo2maxPos;
    std::map<TSeqNo, TSeqPos> seqNo2maxPos;
    for (typename TKLocations::const_iterator it = locations.begin(); it != locations.end(); ++it)
    {
        // resetLimits in mapper may lead to empty kmer occurrences
        if ((it->second).first.size() < FREQ_KMER_MIN)
        {
            std::cout << "continue because FREQ_KMER_MIN = " << FREQ_KMER_MIN << std::endl;
            continue;
        }
        // const auto & [seqNo, seqPos, K] = it->first;
        // if (seqNo2maxPos.find(seqNo) == seqNo2maxPos.end())
        //     seqNo2maxPos[seqNo] = seqPos;
        // else  // update largest kmer position
        //     seqNo2maxPos[seqNo] = std::max(seqNo2maxPos[seqNo], seqPos);

        // it->first contained in list it->second, so we iterate only over list
        TSeqNo seqNo;
        TSeqPos seqPos;
        for (std::vector<TLocation>::const_iterator it_loc_fwd = it->second.first.begin(); it_loc_fwd != it->second.first.end(); ++it_loc_fwd)
        {
            seqNo = seqan::getValueI1<TSeqNo, TSeqPos>(*it_loc_fwd);
            seqPos = seqan::getValueI2<TSeqNo, TSeqPos>(*it_loc_fwd);
            if (seqNo2maxPos.find(seqNo) == seqNo2maxPos.end())
                seqNo2maxPos[seqNo] = seqPos;
            else
                seqNo2maxPos[seqNo] = std::max(seqNo2maxPos[seqNo], seqPos);
        }
    }

    std::set<TSeqNo> seqNos;
    for (auto it = seqNo2maxPos.cbegin(); it != seqNo2maxPos.cend(); ++it)
        seqNos.insert(it->first);

    TSeqNo seqNo_cx = 0;
    // store seqNo -> seqNo_cx and (1<<63 | seqNo_cx) -> seqNo (we never have more than 2^63 sequences)
    for (auto seqNo : seqNos)
    {
        seqNoMap[seqNo] = seqNo_cx++;
        seqNoMap[ONE_LSHIFT_63 | seqNoMap[seqNo]] = seqNo;
        // std::cout << seqNo << " ";
    }
    std::cout << std::endl;

    // (ii) Create bit vectors in the length of largest kmer occurences, and set
    // bits for kmer occurrences. Collect also kmer lengths encoded in final kmer
    // code of type uint64_t. Final sequence lookup and encoding is done in next step.
    // seqNoMap.clear();

    // reserve space for bit vectors
    references.resize(seqNo2maxPos.size());
    // std::cout << "seqNo2maxPos.size() = " << seqNo2maxPos.size() << std::endl;
    for (auto it = seqNo2maxPos.cbegin(); it != seqNo2maxPos.cend(); ++it)
    {
        // std::cout << "resize bit vector for seqNo = " << it->first << std::endl;
        sdsl::bit_vector bv(it->second + 1, 0);
        references[seqNoMap[it->first]] = bv;
    }

    std::cout << "STATUS: bit vector reference assignment done\n";
    // set bits for kmers in references and collect Ks
    // bit 0: set if kmer of length primer_length_min found ...
    // bit primer_length_max-primer_length_min: set if kmer of length primer_length_max found
    // => maximal encodable length difference is 16!
    std::unordered_map<uint64_t, uint64_t> loc2k;
    kmerIDs.clear();
    kmerIDs.resize(references.size());
    std::vector<uint32_t> debug_drop_kmer_repeats(references.size(), 0);
    for (typename TKLocations::const_iterator it = locations.begin(); it != locations.end(); ++it)
    {
        // cleanup in mapper may lead to undercounting kmer occurrences
        // TODO: move kmer frequency cutoff completely into mapper
        if ((it->second).first.size() < FREQ_KMER_MIN)
        {
             std::cout << "WARNING: kmer location sizes undershoots FREQ_KMER_MIN\n";
             continue;
        }
        const auto & [seqNo, seqPos, K] = (it->first);
        // use symmetry and lexicographical ordering of locations to skip already seen ones
        // TODO: is this already shortcut fm mapper?

        if (it->second.first.size() && (seqan::getValueI1<TSeqNo, TSeqPos>(it->second.first[0]) < seqNo ||
            (seqan::getValueI1<TSeqNo, TSeqPos>(it->second.first[0]) == seqNo &&
            seqan::getValueI2<TSeqNo, TSeqPos>(it->second.first[0]) < seqPos)))
        {
            std::cout << "WARNING: symmetry of location ordering not exploited! Fix this in mapper.\n";
            continue;
        }
        // insert now unique occurrences listed in first vector (forward mappings)
        TSeqNo seqNo_prev{0};
        TSeqPos seqPos_prev{0};
        for (std::vector<TLocation>::const_iterator it_loc_fwd = it->second.first.begin(); it_loc_fwd != it->second.first.end(); ++it_loc_fwd)
        {
            // std::cout << "a1\n";
            TSeqNo seqNo = seqan::getValueI1<TSeqNo, TSeqPos>(*it_loc_fwd);
            TSeqPos seqPos = seqan::getValueI2<TSeqNo, TSeqPos>(*it_loc_fwd);
            // if (seqNoMap.find(seqNo) == seqNoMap.end())
            //     std::cout << "ERROR: cannot find seqNo = " << seqNo << " in seqNoMap\n";
            // auto seqNo_cx = seqNoMap.at(seqNo); // compressed sequence id
            // continue and delete previous bit if same kmer occurs within 400 bp
            if (it_loc_fwd > it->second.first.begin() && seqNo_prev == seqNo && seqPos_prev + TRAP_DIST >= seqPos)
            {
                seqan::DnaString seq = seqan::valueById(text, seqNo);
                // undo previous kmer bit in reference
                references[seqNoMap[seqNo]][seqPos_prev] = 0;
                seqPos_prev = seqPos;
                debug_drop_kmer_repeats[seqNo]++;
                continue;
            }
            if (references.size() <= seqNoMap[seqNo] || references[seqNoMap[seqNo]].size() <= seqPos)
            {
                if (references.size() > seqNoMap[seqNo])
                    std::cout << "seqPos = " << seqPos << " and references[seqNo].size() = " << references[seqNoMap[seqNo]].size() << std::endl;
                exit(0);
            }
            references[seqNoMap[seqNo]][seqPos] = 1;
            // std::cout << "a10\n";
            uint64_t loc_key = location_encode(seqNoMap[seqNo], seqPos);
            uint64_t loc_val = ONE_LSHIFT_63 >> (K - PRIMER_MIN_LEN); // later prefix of kmerID
            // locations filled for K_min to K_max
            if (loc2k.find(loc_key) == loc2k.end())
                loc2k.insert({loc_key, loc_val});
            else // update
                loc2k[loc_key] |= loc_val;
            seqNo_prev = seqNo;
            seqPos_prev = seqPos;
        }
    }
    /* print drop out statistics */
    for (unsigned i = 0; i < debug_drop_kmer_repeats.size(); ++i)
    {
        if (debug_drop_kmer_repeats[i])
            std::cerr << "seqNo = " << i << ", dropouts: " << int(debug_drop_kmer_repeats[i]) << std::endl;
    }
    // (iii) lookup kmer sequences, filter and encode as 64 bit integers.
    for (TSeqNo seqNo_cx = 0; seqNo_cx < references.size(); ++seqNo_cx) // Note: we iterate over compressed sequence identifiers
    {
        // std::cout << "STATUS: bit vector assignment for seqNo_cx = " << seqNo_cx << std::endl;
        sdsl::rank_support_v5<1> r1s(&references[seqNo_cx]); // check if once initialized modification of references[i] does not invalidate rank support
        sdsl::select_support_mcl<1,1> s1s(&references[seqNo_cx]);
        for (unsigned r = r1s.rank(references[seqNo_cx].size()); r > 0; --r)
        {
            TSeqPos seqPos = s1s.select(r);
            // TODO: use seqNoMap for all seqNo occurrences
            uint64_t loc_key = location_encode(seqNo_cx, seqPos);

            if (loc2k.find(loc_key) == loc2k.end())
                throw std::invalid_argument("ERROR: " + std::to_string(loc_key) + " not in loc2k dictionary.");

            // get kmerID prefix
            TKmerID kmerID = loc2k[loc_key];
            if (!kmerID)
                throw std::invalid_argument("ERROR: prefix is 0!");

            // identify lowest set bit in prefix
            TKmerLength k_max = PRIMER_MAX_LEN - ffsll(kmerID >> 54) + 1;

            // lookup sequence in corpus and encode
            seqan::DnaString seq = seqan::valueById(text, seqNoMap[ONE_LSHIFT_63 | seqNo_cx]);
            TSeq const & kmer_str = seqan::infixWithLength(seq, seqPos, k_max);
            kmerID |= dna_encoder(kmer_str);

            std::string cs = seqan::toCString(static_cast<seqan::CharString>(kmer_str));

            // erase those length bits in prefix corresponding to kmers not passing the filter
            chemical_filter_single_pass(kmerID);

            // do not store Kmer and reset bit in reference
            if (!(PREFIX_SELECTOR & kmerID))
                references[seqNo_cx][seqPos] = 0;
            else
            {
                kmerCounts->at(KMER_COUNTS::FILTER1_CNT) += __builtin_popcountll(kmerID >> 54);
                kmerIDs[seqNo_cx].push_front(kmerID);
            }
        }
    }
}

/* Combine based on suitable location distances s.t. transcript length is in permitted range.
 * Chemical suitability will be tested by a different function. First position indicates,
 * that the k-mer corresponds to a forward primer, and second position indicates reverse
 * primer, i.e. (k1, k2) != (k2, k1).
 */
template<typename TPairList>
void combine(TReferences const & references, TKmerIDs const & kmerIDs, TPairList & pairs, TKmerCounts * kmerCounts = nullptr)
{
    pairs.clear();
    auto cp_ctr = 0ULL;
    for (uint64_t seqNo_cx = 0; seqNo_cx < references.size(); ++seqNo_cx)
    {
        sdsl::bit_vector reference;
        sdsl::util::assign(reference, references.at(seqNo_cx));
        sdsl::rank_support_v5<1, 1> r1s(&references.at(seqNo_cx));
        sdsl::select_support_mcl<1> s1s(&reference);

        for (uint64_t r_fwd = 1; r_fwd < r1s.rank(reference.size()); ++r_fwd)
        {
            uint64_t idx_fwd = s1s.select(r_fwd);  // text position of r-th k-mer
            TKmerID kmerID_fwd = kmerIDs[seqNo_cx][r_fwd - 1];
            if (!(kmerID_fwd >> CODE_SIZE))
                std::cerr << "ERROR: k length pattern is zero\n";

            // minimal window start position for pairing kmer
            uint64_t w_begin = idx_fwd + PRIMER_MIN_LEN + TRANSCRIPT_MIN_LEN;

            // maximal window end position (exclusive) for pairing kmer
            uint64_t w_end = std::min(reference.size(), idx_fwd + PRIMER_MAX_LEN + TRANSCRIPT_MAX_LEN + 1);

            // iterate through kmers in reference sequence window [w_begin : w_end]
            // note that w_begin/end are updated due to varying kmer length of same kmerID
            for (uint64_t r_rev = r1s.rank(w_begin) + 1; r_rev <= r1s.rank(w_end); ++r_rev)
            {
                TCombinePattern<TKmerID, TKmerLength> cp;
                uint64_t mask_fwd = ONE_LSHIFT_63;
                TKmerID kmerID_rev = kmerIDs.at(seqNo_cx).at(r_rev - 1);
                filter_cross_annealing(kmerID_fwd, kmerID_rev);
                while ((((mask_fwd - 1) << 1) & kmerID_fwd) >> 54)
                {
                    // check forward primer not ending with TTT, ATT
                    if ((mask_fwd & kmerID_fwd) && filter_CG_clamp(kmerID_fwd, '+', mask_fwd) && filter_WWW_tail(kmerID_fwd, '+', mask_fwd))
                    {
                        uint64_t mask_rev = ONE_LSHIFT_63;
                        while ((((mask_rev - 1) << 1) & kmerID_rev) >> 54)
                        {
                            cp_ctr++;
                            if (mask_rev & kmerID_rev && filter_CG_clamp(kmerID_rev, '-') && filter_WWW_tail(kmerID_rev, '-'))
                            {
                                if (dTm(kmerID_fwd, mask_fwd, kmerID_rev, mask_rev) <= PRIMER_DTM)
                                {
                                    // store combination bit
                                    cp.set(mask_fwd, mask_rev);

                                    ++kmerCounts->at(KMER_COUNTS::COMBINER_CNT);
                                }
                            }
                            mask_rev >>= 1; // does not affect search window, since starting position is fixed
                        } // length mask_rev
                    }
                    mask_fwd >>= 1;
                } // length mask_fwd
                if (cp.is_set())
                {
                    pairs.push_back(TPair<TCombinePattern<TKmerID, TKmerLength>>{seqNo_cx, r_fwd, r_rev, cp});
                }
            } // kmerID rev
        } // kmerID fwd
    }
}

// Apply frequency cutoff for unique pair occurences
template<typename TPairLists, typename TResultList>
void filter_and_retransform(TReferences & references, TKmerIDs const & kmerIDs, TPairLists & pairs_grouped, TResultList & results)
{
    std::unordered_map<uint64_t, uint32_t> pairhash2freq;
    std::cout << "STATUS: Enter filter_pairs ...\n";
    // count unique pairs
    for (auto it_pairs = pairs.begin(); it_pairs != pairs.end(); ++it_pairs)
    {
        std::vector<std::pair<uint8_t, uint8_t>> combinations;
        it_pairs->cp.get_combinations(combinations);
        TKmerID kmerID_fwd = kmerIDs.at(it_pairs->reference).at(it_pairs->r_fwd - 1);
        TKmerID kmerID_rev = kmerIDs.at(it_pairs->reference).at(it_pairs->r_rev - 1);

        for (auto comb : combinations)
        {
            auto code_fwd = get_code(kmerID_fwd, ONE_LSHIFT_63 >> comb.first);
            auto code_rev = get_code(kmerID_rev, ONE_LSHIFT_63 >> comb.second);
            auto h = hash_pair(code_fwd, code_rev);
            if (pairhash2freq.find(h) == pairhash2freq.end())
                pairhash2freq[h] = 1;
            else
                pairhash2freq[h]++;
        }
    }

    // Reset bits of kmer prefixes if and in references if prefix is finally empty
    uint64_t ctr_reset = 0;
    std::unordered_set<uint64_t> seen;
    for (auto it_pairs = pairs.begin(); it_pairs != pairs.end(); ++it_pairs)
    {
        TKmerID kmer_fwd = kmerIDs.at(it_pairs->reference).at(it_pairs->r_fwd - 1);
        TKmerID kmer_rev = kmerIDs.at(it_pairs->reference).at(it_pairs->r_rev - 1);
        std::vector<std::pair<uint8_t, uint8_t>> combinations;
        it_pairs->cp.get_combinations(combinations);
        for (auto comb : combinations)
        {
            uint64_t code_fwd = get_code(kmer_fwd, ONE_LSHIFT_63 >> comb.first);
            uint64_t code_rev = get_code(kmer_rev, ONE_LSHIFT_63 >> comb.second);
            uint64_t h = hash_pair(code_fwd, code_rev);
            if (pairhash2freq.find(h) == pairhash2freq.end())
            {
                std::cout << "ERROR: unknown pair hash\n";
                exit(0);
            }
            // delete pair if frequency below FREQ_PAIR_MIN
            if (pairhash2freq.at(h) < FREQ_PAIR_MIN)
            {
                ctr_reset++;
                it_pairs->cp.reset(comb.first, comb.second);
            }
            else
            {
                if (seen.find(h) == seen.end())
                {
                    // TODO: derive container types from template parameter
                    std::tuple<uint64_t, uint64_t, uint64_t, uint64_t> pp{kmer_fwd, ONE_LSHIFT_63 >> comb.first, kmer_rev, ONE_LSHIFT_63 >> comb.second};
                    std::pair<uint32_t, std::tuple<uint64_t, uint64_t, uint64_t, uint64_t>> ppf{pairhash2freq.at(h), pp};
                    pair_freqs.push_back(ppf);
                    seen.insert(h);
                }
                kmerCounts->at(KMER_COUNTS::FILTER2_CNT)++;
            }
        }
    }
    std::cout << "STATUS: Leaving filter_pairs\n";
}

}  // namespace priset
