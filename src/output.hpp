#pragma once

#include <fstream>
#include <iostream>
#include <iterator>
#include <set>
#include <string>
#include <unordered_set>

#include "chemistry.hpp"
#include "types/all.hpp"
#include "utilities.hpp"

namespace priset
{
// forward declaration
uint64_t get_code(uint64_t const code_, uint64_t mask);
std::string dna_decoder(uint64_t code, uint64_t const mask);

//template<typename TKmerIDs>
//void unique_kmers(TKmerIDs const & kmerIDs, std::set<uint64_t> & set);

// Result output helper for writing primer infos.
template<typename IOConfig, typename PrimerConfig, typename TKmerIDs>
void write_primer_info_file(IOConfig const & io_cfg, PrimerConfig const & primer_cfg, TKmerIDs const & kmerIDs)
{
    using TKmerID = typename TKmerIDs::value_type::value_type;
    std::ofstream primer_table;
    primer_table.open(io_cfg.get_primer_info_file());
    primer_table << "id,kmer_sequence,coverage_tax,coverage_ref,length,CG,Tm\n";
    std::set<uint64_t> kmer_ordered_set;
    unique_kmers(kmerIDs, kmer_ordered_set);
    uint64_t i = 0;
    for (TKmerID kmerID : kmer_ordered_set)
    {
        primer_table << (i++) << "," << dna_decoder(kmerID, 0) << "," << get_Tm(primer_cfg, kmerID) << "\n";
    }
    primer_table.close();
    std::cout << "STATUS: primer_info.csv written to\t" << io_cfg.get_primer_info_file() << std::endl;
}

/*
 * Write sequences to file without further sequence information. A handle of the
 * unique kmers allows further processing of the caller.
*/
template<typename TKmerIDs, typename PairList>
void write_primer_file(TKmerIDs const & kmerIDs, PairList const & pairs, fs::path const & primer_file, std::unordered_set<std::string> & kmers_unique_str)
{
    std::unordered_set<uint64_t> kmers_unique;
    for (auto pair : pairs)
    {
        uint64_t kmer_fwd = kmerIDs.at(pair.reference).at(pair.r_fwd - 1);
        uint64_t kmer_rev = kmerIDs.at(pair.reference).at(pair.r_rev - 1);

        for (uint8_t i = 0; i < 100; ++i)
        {
            if (pair.cp[i])
            {
                uint64_t code_fwd = get_code(kmer_fwd & ~PREFIX_SELECTOR, ONE_LSHIFT_63 >> (i/10));
                uint64_t code_rev = get_code(kmer_rev & ~PREFIX_SELECTOR, ONE_LSHIFT_63 >> (i % 10));

                if (!code_fwd)
                {
                    std::cout << "ERROR: kmer_fwd = " << kmer_fwd << " called with get_code(..., " << (ONE_LSHIFT_63 >> (i/10)) << ") is 0ULL\n";
                    exit(0);
                }
                if (!code_rev)
                {
                    std::cout << "ERROR: kmer_fwd = " << kmer_rev << " called with get_code(..., 1<<(63 - " << (i % 10) << ") is 0ULL\n";
                    exit(0);
                }
                kmers_unique.insert(code_fwd);
                kmers_unique.insert(code_rev);
            }
        }
    }

    std::ofstream ofs(std::string(primer_file), std::ofstream::out);
    ofs << "primer\n";
    for (auto primer : kmers_unique)
    {
        auto primer_str = dna_decoder(primer, 0);
        kmers_unique_str.insert(primer_str);
        ofs << primer_str << "\n";
    }
    ofs.close();
    std::cout << "MESSAGE: output written to " << primer_file << std::endl;

}

struct pair_hash
{
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2> &pair) const
    {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};

template<typename PairList, typename TKmerIDs, typename Pair2RefMap>
void unfold_pairs(PairList const & pairs, TKmerIDs const & kmerIDs, Pair2RefMap & pair2refs)
{
    for (auto pair : pairs)
    {
        auto refID = pair.reference;
        auto code1 = ~PREFIX_SELECTOR & kmerIDs.at(refID).at(pair.r_fwd - 1);
        auto code2 = ~PREFIX_SELECTOR & kmerIDs.at(refID).at(pair.r_rev - 1);
        std::vector<std::pair<uint8_t, uint8_t>> combinations;
        pair.cp.get_combinations(combinations);
        for (auto combination : combinations)
        {
            uint64_t code1_trim = get_code(code1, ONE_LSHIFT_63 >> (combination.first));
            uint64_t code2_trim = get_code(code2, ONE_LSHIFT_63 >> (combination.second));
            auto key = std::pair<uint64_t, uint64_t>{code1_trim, code2_trim};
            if (pair2refs.find(key) == pair2refs.end())
                pair2refs[key] = std::vector<uint64_t>{refID};

        }
    }
}

/*
    Write result table with columns: taxid, fwd, rev, matches, coverage, ID_list and
    primer info file with columns kmer_id (1-based), sequence and melting temperature.
    Output clade_x_results.csv:
    fwd         forward primer ID
    rev         reverse primer ID
    dTm         difference in melting temperature
    cov_tax     number of taxa covered (at least one reference per taxon)
    cov_ref     number of references matching
    acc_list    list of accessions that are matched

    Output primers.csv
    id          unique primer ID
    seq         DNA sequence
    len         primer length
    Tm          melting temperature
    CG          CG content as float
*/
template<typename IOConfig, typename PrimerConfig, typename PairList, typename TSeqNoMap, typename TReferences, typename KmerIDs, typename Result>
void create_table(IOConfig const & io_cfg, PrimerConfig const & primer_cfg, /*TSeqNoMap const & seqNoMap,*/TReferences const & references, TKmerIDs const & kmerIDs, PairList const & pairs)
{
    //using TKmerID = typename TKmerIDs::value_type::value_type;
    //std::set<uint64_t> kmer_ordered_set;
    //unique_kmers(kmerIDs, kmer_ordered_set);
    // pair: code fwd, code rev referring to list of reference IDs
    using Pair2RefMap = std::unordered_map<std::pair<uint64_t, uint64_t>, std::vector<uint64_t>, pair_hash>;
    Pair2RefMap pair2refs;
    unfold_pairs<PairList, TKmerIDs, Pair2RefMap>(pairs, kmerIDs, pair2refs);

    //std::cout << "STATUS: will write " << kmer_ordered_set.size() << " single primer results and " << pairs.size() << " primer pair results\n";
    std::ofstream table;
    table.open(io_cfg.get_result_file());
    std::cout << "write header line\n";
    table << "fwd,rev,dTm,cov_tax,cov_ref,acc_list\n"; // no taxid here!
    table.close();

/*
    // TODO: check if unordered_map instead of map
    // load id file for mapping reference IDs (1-based) to accession numbers and vice versa
    std::unordered_map<AccessionID, Accession> accID2acc;
    std::unordered_map<Accession, AccessionID> acc2accID;
    create_accID2acc_map<IOConfig>(accID2acc, acc2accID, io_cfg);

    // build dictionary for taxids and counter for assigned accessions
    std::unordered_map<AccessionID, Taxid> accID2taxID;
    std::unordered_set<Taxid> taxid_set; // taxids with accessions
    create_accID2taxID_map(accID2taxID, taxid_set, acc2accID, io_cfg);

    // load taxonomy as map {taxid: p_taxid}, taxid is root if not in key set
    std::unordered_map<Taxid, Taxid> tax_map;
    create_tax_map(tax_map, io_cfg);

    // collect single kmer matches for bottom nodes
    std::unordered_map<TKmerID, std::vector<TSeqNo> > kmer2loc; // relates kmer IDs and location IDs
    for (auto it = kmer_locations.begin(); it != kmer_locations.end(); ++it)
    {
        TKmerID kmer_ID = it->get_kmer_ID();
        std::vector<TSeqNo> seq_IDs;
        for (KmerLocation::size_type i = 0; i < it->container_size(); ++i)
            seq_IDs.push_back(it->accession_ID_at(i)); // seqan::getValueI1<TSeqNo, TSeqPos>(loc));
        kmer2loc[kmer_ID] = seq_IDs;
    }
    // 0-based height, correct level info iff a taxonomic node is in the predecessor lineage of another one
    std::unordered_map<Taxid, uint16_t> leaves; //<taxid, height_from_bottom>

    for (auto const & taxid : taxid_set)
        leaves[taxid] = 0;
    for (auto leaf_it = leaves.begin(); leaf_it != leaves.end(); ++leaf_it)
    {
        uint16_t level = leaf_it->second;
        Taxid p_taxid = leaf_it->first;
        while (tax_map.find(p_taxid) != tax_map.end()) //auto tax_it{tax_map.find(p_taxid)} != tax_map.end())
        {
            ++level;
            p_taxid = tax_map[p_taxid]; //tax_it.second; // set taxid to ancestor node

            auto anc_it = leaves.find(p_taxid);
            if (anc_it != leaves.end()) // && !(p_taxid < (*leaf_it).first)))
            {
                anc_it->second = std::max<uint16_t>((anc_it->second), level);
            }
        }
    }

    // taxids sorted by level to have correct upstreams stats
    std::vector<std::pair<Taxid, uint16_t>> leaves_srt_by_level(leaves.size());
    std::copy(leaves.begin(), leaves.end(), leaves_srt_by_level.begin());
    std::sort(leaves_srt_by_level.begin(), leaves_srt_by_level.end(), [](auto const & l1, auto const & l2){ return l1.second < l2.second; });

    // write result table header
    std::ofstream table;
    table.open(io_cfg.get_result_file());
    // taxid, fwd primer ID, rev primer ID, number of matches, coverage (ctr of nodes with accessions), accession list (comma separated string)
    std::cout << "write header line\n";
    table << "taxid,fwd,rev,matches,coverage,accession_list\n";
    table.close();

    // collect single kmer matches
    accumulation_loop<KmerLocations, IOConfig, Result>(kmer_locations, leaves_srt_by_level, tax_map, accID2taxID, accID2acc, io_cfg);

    // collect kmer pair matches for bottom nodes
    accumulation_loop<Pairs>(kmer_pairs, leaves_srt_by_level, tax_map, accID2taxID, accID2acc, io_cfg);
    std::cout << "STATUS: table.csv written to\t" << io_cfg.get_result_file() << std::endl;

    // write primer info file
    write_primer_info_file(io_cfg, primer_cfg, kmer_locations);
    */
}

} // namespace priset
