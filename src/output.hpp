#pragma once

#include <fstream>
#include <iostream>
#include <iterator>
#include <set>
#include <string>

#include "chemistry.hpp"
#include "combine_types.hpp"
#include "io_cfg_type.hpp"
#include "primer_cfg_type.hpp"
#include "types.hpp"
#include "utilities.hpp"

// forward declaration
//template<typename TKmerIDs>
//void unique_kmers(TKmerIDs const & kmerIDs, std::set<uint64_t> & set);

// Result output helper for writing primer infos.
template<typename io_cfg_type, typename primer_cfg_type, typename TKmerIDs>
void write_primer_info_file(io_cfg_type const & io_cfg, primer_cfg_type const & primer_cfg, TKmerIDs const & kmerIDs)
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
        primer_table << (i++) << "," << dna_decoder(kmerID) << "," << get_Tm(primer_cfg, kmerID) << "\n";
    }
    primer_table.close();
    std::cout << "STATUS: primer_info.csv written to\t" << io_cfg.get_primer_info_file() << std::endl;
}

struct pair_hash
{
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2> &pair) const
    {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};

template<typename TPairList, typename TKmerIDs, typename TPair2RefMap>
void unfold_pairs(TPairList const & pairs, TKmerIDs const & kmerIDs, TPair2RefMap & pair2refs)
{
    for (auto pair : pairs)
    {
        auto refID = pair.reference;
        auto code1 = ~PREFIX_SELECTOR & kmerIDs.at(refID).at(pair.r_fwd);
        auto code2 = ~PREFIX_SELECTOR & kmerIDs.at(refID).at(pair.r_rev);

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
template<typename io_cfg_type, typename primer_cfg_type, typename TPairList, typename TSeqNoMap, typename TReferences, typename TKmerIDs>
void create_table(io_cfg_type const & io_cfg, primer_cfg_type const & primer_cfg, TSeqNoMap const & seqNoMap, TReferences const & references, TKmerIDs const & kmerIDs, TPairList const & pairs)
{
    //using TKmerID = typename TKmerIDs::value_type::value_type;
    //std::set<uint64_t> kmer_ordered_set;
    //unique_kmers(kmerIDs, kmer_ordered_set);
    // pair: code fwd, code rev referring to list of reference IDs
    using TPair2RefMap = std::unordered_map<std::pair<uint64_t, uint64_t>, std::vector<uint64_t>, pair_hash>;
    TPair2RefMap pair2refs;
    unfold_pairs<TPairList, TKmerIDs, TPair2RefMap>(pairs, kmerIDs, pair2refs);

    //std::cout << "STATUS: will write " << kmer_ordered_set.size() << " single primer results and " << pairs.size() << " primer pair results\n";
    std::ofstream table;
    table.open(io_cfg.get_result_file());
    std::cout << "write header line\n";
    table << "fwd,rev,dTm,cov_tax,cov_ref,acc_list\n"; // no taxid here!
    table.close();

/*
    // TODO: check if unordered_map instead of map
    // load id file for mapping reference IDs (1-based) to accession numbers and vice versa
    std::unordered_map<TAccID, TAcc> accID2acc;
    std::unordered_map<TAcc, TAccID> acc2accID;
    create_accID2acc_map<io_cfg_type>(accID2acc, acc2accID, io_cfg);

    // build dictionary for taxids and counter for assigned accessions
    std::unordered_map<TAccID, TTaxid> accID2taxID;
    std::unordered_set<TTaxid> taxid_set; // taxids with accessions
    create_accID2taxID_map(accID2taxID, taxid_set, acc2accID, io_cfg);

    // load taxonomy as map {taxid: p_taxid}, taxid is root if not in key set
    std::unordered_map<TTaxid, TTaxid> tax_map;
    create_tax_map(tax_map, io_cfg);

    // collect single kmer matches for bottom nodes
    std::unordered_map<TKmerID, std::vector<TSeqNo> > kmer2loc; // relates kmer IDs and location IDs
    for (auto it = kmer_locations.begin(); it != kmer_locations.end(); ++it)
    {
        TKmerID kmer_ID = it->get_kmer_ID();
        std::vector<TSeqNo> seq_IDs;
        for (TKmerLocation::size_type i = 0; i < it->container_size(); ++i)
            seq_IDs.push_back(it->accession_ID_at(i)); // seqan::getValueI1<TSeqNo, TSeqPos>(loc));
        kmer2loc[kmer_ID] = seq_IDs;
    }
    // 0-based height, correct level info iff a taxonomic node is in the predecessor lineage of another one
    std::unordered_map<TTaxid, uint16_t> leaves; //<taxid, height_from_bottom>

    for (auto const & taxid : taxid_set)
        leaves[taxid] = 0;
    for (auto leaf_it = leaves.begin(); leaf_it != leaves.end(); ++leaf_it)
    {
        uint16_t level = leaf_it->second;
        TTaxid p_taxid = leaf_it->first;
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
    std::vector<std::pair<TTaxid, uint16_t>> leaves_srt_by_level(leaves.size());
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
    accumulation_loop<TKmerLocations, io_cfg_type>(kmer_locations, leaves_srt_by_level, tax_map, accID2taxID, accID2acc, io_cfg);

    // collect kmer pair matches for bottom nodes
    accumulation_loop<TPairs>(kmer_pairs, leaves_srt_by_level, tax_map, accID2taxID, accID2acc, io_cfg);
    std::cout << "STATUS: table.csv written to\t" << io_cfg.get_result_file() << std::endl;

    // write primer info file
    write_primer_info_file(io_cfg, primer_cfg, kmer_locations);
    */
}
