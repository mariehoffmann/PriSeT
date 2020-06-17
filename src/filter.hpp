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
// #include "../submodules/sdsl-lite/include/sdsl/bit_vectors.hpp"

#include "types/all.hpp"
#include "utilities.hpp"

namespace priset
{

// Filter of single kmers and transform of references to bit vectors.
template<typename TKLocations, typename TSeqNoMap, typename TKmerIDs>
void transform_and_filter(IOConfig const & io_cfg, PrimerConfig & primer_cfg, TKLocations & locations, TReferences & references, TSeqNoMap & seqNo_map, TKmerIDs & kmerIDs, std::vector<Taxid> & taxa_by_seqNo_cx, uint64_t * kmerCounts = nullptr)
{
    // uniqueness indirectly preserved by (SeqNo, SeqPos) if list sorted lexicographically
    assert(length(locations));

    references.clear();
    kmerIDs.clear();
    taxa_by_seqNo_cx.clear();

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
        if ((it->second).first.size() < primer_cfg.get_digamma())
        {
            std::cout << "continue because digamma = " << primer_cfg.get_digamma() << std::endl;
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
    for (TSeqNo seqNo : seqNos)
    {
        seqNo_map[seqNo] = seqNo_cx++;
        seqNo_map[ONE_LSHIFT_63 | seqNo_map[seqNo]] = seqNo;
        // std::cout << seqNo << " ";
        // TODO: is seqNo 0-based?
        Accession acc = io_cfg.get_seqNo_by_accession(seqNo); //seqNo2acc_map.at(seqNo);
        // [taxid_seq0, taxid_seq1, ...]
        taxa_by_seqNo_cx.push_back(io_cfg.get_taxid_by_accession(acc)); //acc2taxid_map.at(acc));
    }
    std::cout << std::endl;

    // (ii) Create bit vectors in the length of largest kmer occurences, and set
    // bits for kmer occurrences. Collect also kmer lengths encoded in final kmer
    // code of type uint64_t. Final sequence lookup and encoding is done in next step.
    // seqNo_map.clear();

    // reserve space for bit vectors
    references.resize(seqNo2maxPos.size());
    // std::cout << "seqNo2maxPos.size() = " << seqNo2maxPos.size() << std::endl;
    for (auto it = seqNo2maxPos.cbegin(); it != seqNo2maxPos.cend(); ++it)
    {
        // std::cout << "resize bit vector for seqNo = " << it->first << std::endl;
        sdsl::bit_vector bv(it->second + 1, 0);
        references[seqNo_map[it->first]] = bv;
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
        if ((it->second).first.size() < primer_cfg.get_digamma())
        {
             std::cout << "WARNING: kmer location sizes undershoots digamma\n";
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
            // if (seqNo_map.find(seqNo) == seqNo_map.end())
            //     std::cout << "ERROR: cannot find seqNo = " << seqNo << " in seqNo_map\n";
            // auto seqNo_cx = seqNo_map.at(seqNo); // compressed sequence id
            // continue and delete previous bit if same kmer occurs within 400 bp
            if (it_loc_fwd > it->second.first.begin() && seqNo_prev == seqNo && seqPos_prev + TRAP_DIST >= seqPos)
            {
                seqan::DnaString seq = seqan::valueById(text, seqNo);
                // undo previous kmer bit in reference
                references[seqNo_map[seqNo]][seqPos_prev] = 0;
                seqPos_prev = seqPos;
                debug_drop_kmer_repeats[seqNo]++;
                continue;
            }
            if (references.size() <= seqNo_map[seqNo] || references[seqNo_map[seqNo]].size() <= seqPos)
            {
                if (references.size() > seqNo_map[seqNo])
                    std::cout << "seqPos = " << seqPos << " and references[seqNo].size() = " << references[seqNo_map[seqNo]].size() << std::endl;
                exit(0);
            }
            references[seqNo_map[seqNo]][seqPos] = 1;
            // std::cout << "a10\n";
            uint64_t loc_key = location_encode(seqNo_map[seqNo], seqPos);
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
            // TODO: use seqNo_map for all seqNo occurrences
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
            seqan::DnaString seq = seqan::valueById(text, seqNo_map[ONE_LSHIFT_63 | seqNo_cx]);
            TSeq const & kmer_str = seqan::infixWithLength(seq, seqPos, k_max);
            kmerID |= dna_encoder(kmer_str);

            std::string cs = seqan::toCString(static_cast<seqan::CharString>(kmer_str));

            // erase those length bits in prefix corresponding to kmers not passing the filter
            chemical_filter_single_pass(kmerID);

            // do not store Kmer and reset bit in reference
            if (!(PREFIX_SELECTOR & kmerID))
                references[seqNo_cx][seqPos] = 0;
            else
                kmerIDs[seqNo_cx].push_front(kmerID);
        }
    }
}

/* Combine based on suitable location distances s.t. transcript length is in permitted range.
 * Chemical suitability will be tested by a different function. First position indicates,
 * that the k-mer corresponds to a forward primer, and second position indicates reverse
 * primer, i.e. (k1, k2) != (k2, k1).
 */
template<typename PrimerPairList, typename TKmerIDs>
void combine(TReferences const & references, TKmerIDs const & kmerIDs, PrimerPairList & pairs, uint64_t * kmer_counts = nullptr)
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
                CombinePattern cp;
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
                                    ++kmer_counts[STEPS::COMBINE];
                                }
                            }
                            mask_rev >>= 1; // does not affect search window, since starting position is fixed
                        } // length mask_rev
                    }
                    mask_fwd >>= 1;
                } // length mask_fwd
                if (cp.is_set())
                {
                    pairs.push_back(PrimerPair{seqNo_cx, r_fwd, r_rev, cp});
                }
            } // kmerID rev
        } // kmerID fwd
    }
}

// Filter pairs by frequency and unpack PrimerPair -> PrimerPairUnpacked
template<typename TSeqNoMap, typename TKmerIDs, typename PrimerPairList, typename PrimerPairUnpackedList>
void filter_and_unpack_pairs(IOConfig & io_cfg, PrimerConfig & primer_cfg, TSeqNoMap & seqNo_map, TKmerIDs const & kmerIDs, PrimerPairList & pairs, PrimerPairUnpackedList & pairs_unpacked)
{
    std::unordered_map<std::string, std::vector<bool>> pair2seqNo_cx_vector;
    std::unordered_map<uint64_t, std::string> primers_memoized;
    std::vector<std::pair<uint8_t, uint8_t>> combinations;

    auto pair_key = [](std::string const & p1, std::string const & p2) -> std::string
        {return p1 + "_" + p2;};
    std::cout << "STATUS: Enter filter_pairs ...\n";
    // count unique pairs
    for (auto it_pairs = pairs.begin(); it_pairs != pairs.end(); ++it_pairs)
    {
        TSeqNo seqNo_cx = it_pairs->get_seqNo();
        it_pairs->get_combine_pattern().get_combinations(combinations);
        TKmerID kmerID_fwd = kmerIDs.at(it_pairs->get_seqNo()).at(it_pairs->get_rank_fwd() - 1);
        TKmerID kmerID_rev = kmerIDs.at(it_pairs->get_seqNo()).at(it_pairs->get_rank_rev() - 1);

        for (auto comb : combinations)
        {
            uint64_t code_fwd = get_code(kmerID_fwd, ONE_LSHIFT_63 >> comb.first);
            uint64_t code_rev = get_code(kmerID_rev, ONE_LSHIFT_63 >> comb.second);
            std::string primer_fwd = dna_decoder(kmerID_fwd, code_fwd);
            std::string primer_rev = dna_decoder(kmerID_rev, code_rev);

            if (!primers_memoized.count(code_fwd))
                primers_memoized[code_fwd] = primer_fwd;
            if (!primers_memoized.count(code_rev))
                primers_memoized[code_rev] = primer_rev;

            std::string key = pair_key(primer_fwd, primer_rev);
            if (!pair2seqNo_cx_vector.count(key))
                pair2seqNo_cx_vector[key] = std::vector<bool>(0, seqNo_cx);
            else if (seqNo_cx >= pair2seqNo_cx_vector.at(key).size())
                pair2seqNo_cx_vector[key].resize(0, seqNo_cx + 1);
            pair2seqNo_cx_vector[key][seqNo_cx] = 1;
        }
    }

    // Unpack pairs exceeding minimal pair frequency (digamma_pairs) and store in
    // pairs_unpacked. Combination bits are not reset.
    std::unordered_map<std::string, size_t> seen_and_index;
    for (auto it_pairs = pairs.begin(); it_pairs != pairs.end(); ++it_pairs)
    {
        TKmerID kmer_fwd = kmerIDs.at(it_pairs->get_seqNo()).at(it_pairs->get_rank_fwd() - 1);
        TKmerID kmer_rev = kmerIDs.at(it_pairs->get_seqNo()).at(it_pairs->get_rank_rev() - 1);
        it_pairs->get_combine_pattern().get_combinations(combinations);
        for (auto comb : combinations)
        {
            uint64_t code_fwd = get_code(kmer_fwd, ONE_LSHIFT_63 >> comb.first);
            uint64_t code_rev = get_code(kmer_rev, ONE_LSHIFT_63 >> comb.second);

            std::string key = pair_key(primers_memoized[code_fwd], primers_memoized[code_rev]);
            size_t bit_count = std::accumulate(pair2seqNo_cx_vector.at(key).begin(), pair2seqNo_cx_vector.at(key).end(), 0);
            if (bit_count >= primer_cfg.get_digamma_pairs())
            {
                if (!seen_and_index.count(key))
                {
                    PrimerPairUnpacked<TSeqNoMap> pair_unpacked{&io_cfg, &seqNo_map, pair2seqNo_cx_vector.at(key), kmer_fwd, code_fwd, kmer_rev, code_rev};
                    pairs_unpacked.push_back({pair_unpacked});
                    seen_and_index[key] = pairs_unpacked.size() - 1;
                }
                pairs_unpacked.at(seen_and_index.at(key)).set_sequence_match(it_pairs->get_seqNo());
            }
        }
    }
}

// Apply frequency cutoff for unique pair occurences and unfold KmerIDs to DNA sequences.
// template<typename TKmerIDs, typename PrimerPairGroups, typename ResultList>
// void filter_groups(PrimerConfig const & primer_cfg, TReferences & references, TKmerIDs const & kmerIDs, PrimerPairGroups & pairs_grouped, ResultList & results, uint64_t * kmerCounts = nullptr)
// {
//
//     for (PrimerPairGroups::value_type group : pairs_grouped)
//     {
//
//
//
//
//     }
//     std::cout << "STATUS: Leaving filter_pairs\n";
// }

}  // namespace priset
