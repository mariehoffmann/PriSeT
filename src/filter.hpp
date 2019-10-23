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

// Filter of single kmers and transform of references to bit vectors.
void filter_and_transform(io_cfg_type const & io_cfg, TKLocations & locations, TReferences & references, TKmerIDs & kmerIDs, TSeqNoMap & seqNoMap, TKmerCounts & stats)
{
    // uniqueness indirectly preserved by (SeqNo, SeqPos) if list sorted lexicographically
    assert(length(locations));

    references.clear();
    kmerIDs.clear();
    // lower bound for kmer frequency
    unsigned freq_kmer_min_abs = FREQ_KMER_MIN; //io_cfg.get_freq_kmer_min();

    // load corpus for dna to 64 bit conversion
    seqan::StringSet<seqan::DnaString, seqan::Owner<seqan::ConcatDirect<>>> text;
    fs::path text_path = io_cfg.get_index_txt_path();
    seqan::open(text, text_path.string().c_str(), seqan::OPEN_RDONLY);
    if (!seqan::lengthSum(text))
        throw std::length_error("Reference text size is 0!");
    // (i) collect distinct sequence identifiers and maximal position of kmer occurences
    // to have a compressed representation.
    std::map<TSeqNo, TSeqPos> seqNo2maxPos;
    for (typename TKLocations::const_iterator it = locations.begin(); it != locations.end(); ++it)
    {
        // resetLimits in mapper may lead to empty kmer occurrences
        if ((it->second).first.size() < freq_kmer_min_abs)
        {
            std::cout << "continue because FREQ_KMER_MIN = " << freq_kmer_min_abs << std::endl;
            continue;
        }
        const auto & [seqNo, seqPos, K] = it->first;
        // use symmetry and lexicographical ordering of locations to skip already seen ones
        // TODO: is this already shortcut fm mapper? TKLoc1 <= TLoc2?
        if (it->second.first.size() && (std::get<0>(it->first) > it->second.first[0].i1 ||
           (std::get<0>(it->first) == it->second.first[0].i1 && std::get<1>(it->first) > it->second.first[0].i2)))
        {
            std::cout << "WARNING: symmetry of location ordering not exploited! Fix this in mapper.\n";
            std::cout << "current seqNo = " << seqNo << ", seqPos = " << seqPos << std::endl;
            std::cout << "first entry: seqNo[0] = " << seqan::getValueI1<TSeqNo, TSeqPos>(it->second.first[0]) << ", seqPos[0] = " << seqan::getValueI2<TSeqNo, TSeqPos>(it->second.first[0]) << std::endl;
            continue;
        }
        // update largest kmer position
        seqNo2maxPos[seqNo] = (seqNo2maxPos.find(seqNo) == seqNo2maxPos.end()) ? seqPos : std::max(seqNo2maxPos.at(seqNo), seqPos);
        // and also for other occurences
        for (std::vector<TLocation>::const_iterator it_loc_fwd = it->second.first.begin(); it_loc_fwd != it->second.first.end(); ++it_loc_fwd)
        {
            TSeqNo seqNo_fwd = seqan::getValueI1<TSeqNo, TSeqPos>(*it_loc_fwd);
            TSeqPos seqPos_fwd = seqan::getValueI2<TSeqNo, TSeqPos>(*it_loc_fwd);
            seqNo2maxPos[seqNo_fwd] = (seqNo2maxPos.find(seqNo_fwd) == seqNo2maxPos.end()) ? seqPos_fwd : std::max(seqNo2maxPos.at(seqNo_fwd), seqPos_fwd);
        }
    }
    // (ii) Create bit vectors in the length of largest kmer occurences, and set
    // bits for kmer occurrences. Collect also kmer lengths encoded in final kmer
    // code of type uint64_t. Final sequence lookup and encoding is done in next step.
    seqNoMap.clear();

    unsigned i_cx = 0;
    // reserve space for bit vectors
    for (auto it = seqNo2maxPos.begin(); it != seqNo2maxPos.end(); ++it)
    {
        seqNoMap[it->first] = i_cx++;  // assign reference index to its compressed index
        sdsl::bit_vector bv(it->second + 1, 0);
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
        //std::cout << "1\n";
        // cleanup in mapper may lead to undercounting kmer occurrences
        // TODO: move kmer frequency cutoff completely into mapper
        if ((it->second).first.size() < freq_kmer_min_abs)
            continue; //std::cout << "WARNING: kmer location sizes undershoots FREQ_KMER_MIN\n";
        const auto & [seqNo, seqPos, K] = (it->first);
        // use symmetry and lexicographical ordering of locations to skip already seen ones
        // TODO: is this already shortcut fm mapper?
        //std::cout << "2\n";
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
        //    std::cout << "a\n";
            TSeqNo seqNo = seqan::getValueI1<TSeqNo, TSeqPos>(*it_loc_fwd);
            TSeqPos seqPos = seqan::getValueI2<TSeqNo, TSeqPos>(*it_loc_fwd);
            if (seqNoMap.find(seqNo) == seqNoMap.end())
                std::cout << "ERROR: cannot find seqNo = " << seqNo << " in seqNoMap\n";
            auto seqNo_cx = seqNoMap.at(seqNo); // compressed sequence id
            // continue and delete previous bit if same kmer occurs within 400 bp
            if (it_loc_fwd > it->second.first.begin() && seqNo_prev == seqNo && seqPos_prev + TRAP_DIST < seqPos)
            {
                seqan::DnaString seq = seqan::valueById(text, seqNo);
                if (seqan::infixWithLength(seq, seqPos, 24) == "GTAGGTGAACCTGCAGAAGGATCA")
                    std::cout << "DEBUG: TRAP_DIST! for V9 rev\n";

                // undo previous kmer bit in reference
                references[seqNo_cx][seqPos_prev] = 0;
                seqPos_prev = seqPos;
                debug_drop_kmer_repeats[seqNo_cx]++;
                continue;
            }
            references[seqNo_cx][seqPos] = 1;
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

    // (iii) lookup kmer sequences, filter and encode as 64 bit integers.
    // inverse reference ID map
    // todo: use one bit vector instead with rank support
    std::unordered_map<TSeqNo, TSeqNo> seqNoMap_inv;
    for (auto it = seqNoMap.begin(); it != seqNoMap.end(); ++it)
        seqNoMap_inv.insert({it->second, it->first});
    kmerIDs.resize(references.size());

    for (unsigned i = 0; i < references.size(); ++i)
    {
        sdsl::rank_support_v5<1> r1s(&references[i]); // check if once initialized modification of references[i] does not invalidate rank support
        sdsl::select_support_mcl<1,1> s1s(&references[i]);
        for (unsigned r = r1s.rank(references[i].size()); r > 0; --r)
        {
            TSeqNo seqNo = seqNoMap_inv[i];
            TSeqPos seqPos = s1s.select(r);
            uint64_t loc_key = location_encode(seqNo, seqPos);

            if (loc_and_ks.find(loc_key) == loc_and_ks.end())
                throw std::invalid_argument("ERROR: " + std::to_string(loc_key) + " not in loc_and_ks dictionary.");

            TKmerID kmerID_prefix = loc_and_ks[loc_key]; // tail not filled yet & CODE_SIZEgth_mask;
            if (!kmerID_prefix)
                throw std::invalid_argument("ERROR: prefix is 0!");
            // identify lowest set bit in head
            TKmerLength k_max = PRIMER_MAX_LEN - ffsll(kmerID_prefix >> 54) + 1;

            // lookup sequence in corpus and encode
            seqan::DnaString seq = seqan::valueById(text, seqNo);

            // Check if bit is set for V9 forward: yes -> does it remain set after filter
            // No -> why not? all length bits until k_max should be set before filtering!
            TSeq const & kmer_str = seqan::infixWithLength(seq, seqPos, k_max);
            TKmerID kmerID = kmerID_prefix + dna_encoder(kmer_str);

            // erase those length bit in prefix corresponding to kmers not passing the filter
            // TODO: kmerID trimmed to maximal encoded length?
            chemical_filter_single_pass(kmerID);

            // do not store Kmer and reset bit in reference
            if (!(PREFIX_SELECTOR & kmerID))
                references[i][seqPos] = 0;
            else
            {
                stats[KMER_COUNTS::FILTER1_CNT] += __builtin_popcountll(kmerID >> 54);
                kmerIDs[i].push_front(kmerID);
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
void combine(TReferences const & references, TKmerIDs const & kmerIDs, TPairList & pairs, TKmerCounts & stats)
{
    pairs.clear();
    auto cp_ctr = 0ULL;
    for (uint64_t i = 0; i < references.size(); ++i)
    {
        sdsl::bit_vector reference;
        sdsl::util::assign(reference, references.at(i));
        sdsl::rank_support_v5<1, 1> r1s(&references.at(i)); // replace after bugfix with
        sdsl::select_support_mcl<1> s1s(&reference);

        //std::cout << "DEBUG: iterate over " << r1s.rank(reference.size()) << " kmer IDs\n";
        for (uint64_t r_fwd = 1; r_fwd < r1s.rank(reference.size()); ++r_fwd)
        {
            uint64_t idx_fwd = s1s.select(r_fwd);  // text position of r-th k-mer
            TKmerID const kmerID_fwd = kmerIDs[i][r_fwd - 1];
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
                while ((((mask_fwd - 1) << 1) & kmerID_fwd) >> 54)
                {
                    // check forward primer not ending with TTT, ATT
                    if ((mask_fwd & kmerID_fwd) && filter_CG_clamp(kmerID_fwd, '+', mask_fwd) && filter_WTT_tail(kmerID_fwd, mask_fwd, '+'))
                    {
                        TKmerID const kmerID_rev = kmerIDs.at(i).at(r_rev - 1);
                        uint64_t mask_rev = ONE_LSHIFT_63;

                        while ((((mask_rev - 1) << 1) & kmerID_rev) >> 54)
                        {
                            cp_ctr++;
                            if (cp_ctr == 1ULL << 63)
                            {
                                std::cout << "cp_ctr reached 1 << 63, exit\n";
                                exit(0);
                            }
                            if (mask_rev & kmerID_rev && filter_CG_clamp(kmerID_rev, '-') && filter_WTT_tail(kmerID_rev, '-'))
                            {
                                if (dTm(kmerID_fwd, mask_fwd, kmerID_rev, mask_rev) <= PRIMER_DTM)
                                {
                                    // store combination bit
                                    cp.set(mask_fwd, mask_rev);

                                    ++stats[KMER_COUNTS::COMBINER_CNT];
                                }
                            }
                            mask_rev >>= 1; // does not affect search window, since starting position is fixed
                        } // length mask_rev
                    }
                    mask_fwd >>= 1;
                } // length mask_fwd
                if (cp.is_set())
                {
                    pairs.push_back(TPair<TCombinePattern<TKmerID, TKmerLength>>{i, r_fwd, r_rev, cp});

                }
            } // kmerID rev
        } // kmerID fwd
    }
    std::cout << "leaving combine\n";
}

// Apply frequency cutoff for unique pair occurences
template<typename TPairList, typename TPairFreq>
void filter_pairs(TReferences & references, TKmerIDs const & kmerIDs, TPairList & pairs, std::vector<TPairFreq> & pair_freqs, TKmerCounts & kmerCounts)
{
    std::unordered_map<uint64_t, uint32_t> pairhash2freq;
    std::cout << "Enter filter_pairs ...\n";
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
                    std::tuple<uint64_t, uint64_t, uint64_t, uint64_t> pp{kmer_fwd, ONE_LSHIFT_63 >> comb.first,
                                                            kmer_rev, ONE_LSHIFT_63 >> comb.second};
                    std::pair<uint32_t, std::tuple<uint64_t, uint64_t, uint64_t, uint64_t>> ppf{pairhash2freq.at(h), pp};
                    pair_freqs.push_back(ppf);
                    seen.insert(h);
                }
                kmerCounts[KMER_COUNTS::FILTER2_CNT]++;
            }
        }
    }
    std::cout << "Leaving filter_pairs\n";
}

}  // namespace priset
