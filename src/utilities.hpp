#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include <seqan/basic.h>

#include "../submodules/genmap/src/genmap_helper.hpp"

#include "chemistry.hpp"
#include "io_cfg_type.hpp"
#include "types.hpp"

namespace priset
{

void print_locations(TLocations & locations)
{
    using key_type = typename TLocations::key_type;
    using TSeqNo = typename seqan::Value<key_type, 1>::Type;
    using TSeqPos = typename seqan::Value<key_type, 2>::Type;
    //using value1_type = typename seqan::Value<typename TLocations::value_type, 1>::Type;
    std::cout << "(SeqNo, SeqPos): [(SeqNo, SeqPos)], [(SeqNo, SeqPos)]\n";
    for (typename TLocations::const_iterator it = locations.begin(); it != locations.end(); ++it)
    {
        std::cout << "(" << seqan::getValueI1<TSeqNo, TSeqPos>(it->first) << ", " << seqan::getValueI2<TSeqNo, TSeqPos>(it->first) << "): [";
        for (seqan::Pair<TSeqNo, TSeqPos> pair : (it->second).first)
            std::cout << "(" << seqan::getValueI1<TSeqNo, TSeqPos>(pair) << ", " << seqan::getValueI2<TSeqNo, TSeqPos>(pair) << ") ";
        std::cout << "], [";
        for (seqan::Pair<TSeqNo, TSeqPos> pair : (it->second).second)
            std::cout << "(" << seqan::getValueI1(pair) << ", " << seqan::getValueI2(pair) << ") ";
        std::cout << "]\n";
    }
}

void print_kmer_locations(TKmerLocations const & kmer_locations, TKmerMap const & kmer_map)
{
    for (typename TKmerLocations::const_iterator it = kmer_locations.begin(); it != kmer_locations.end(); ++it)
    {
        const TKmer kmer{kmer_map.at((*it).first)};
        std::cout << kmer.seq << ": [";
        for (seqan::Pair<priset::TSeqNo, priset::TSeqPos> loc : it->second)
            std::cout << "(" << seqan::getValueI1<TSeqNo, TSeqPos>(loc) << ", " << seqan::getValueI2<TSeqNo, TSeqPos>(loc) << ") ";

        std::cout << "]" << std::endl;
    }
}

void print_pairs(TKmerPairs const & kmer_pairs, TKmerMap const & kmer_map)
{
    std::set<TKmerID> legend;
    std::cout << "\n(kmer ID1, kmer ID2) | reference ID | (pos1, pos2)\n";
    std::cout << "-------------------------------------------------------\n";
    if (!kmer_pairs.size())
        std::cout << "<None>\n";
    for (typename TKmerPairs::const_iterator it = kmer_pairs.begin(); it != kmer_pairs.end(); ++it)
    {
        std::cout << "(" << (*it).kmer_fwd << ", " << (*it).kmer_rev << ")\t\t| ";
        legend.insert((*it).kmer_fwd);
        legend.insert((*it).kmer_rev);
        for (auto it_loc = (*it).pair_locations.begin(); it_loc != (*it).pair_locations.end(); ++it_loc)
            std::cout << std::get<0>(*it_loc) << "\t\t| (" << std::get<1>(*it_loc) << ", " << std::get<2>(*it_loc) << ")\n";
    }
    if (kmer_pairs.size())
    {
        std::cout << "\n\nkmer ID\t | sequence\n------------------------------\n";
        for (TKmerID kmer_ID : legend)
            std::cout << kmer_ID << "\t| " << kmer_map.at(kmer_ID).seq << std::endl;
    }
}

// forward declaration
struct primer_cfg_type;

// Retrieve DNA sequence from txt.concat given a set of locations
// lookup_sequences<primer_cfg_type>(kmer_locations, io_cfg, primer_cfg, directoryInformation);
// locations: [(ID, kmer_locations)]
template<typename primer_cfg_type>
void lookup_sequences(TKmerLocations & kmer_locations, TKmerMap & kmer_map, io_cfg_type const & io_cfg, primer_cfg_type const & primer_cfg, TDirectoryInformation const & directoryInformation)
{
    std::string fastaFile = std::get<0>(retrieveDirectoryInformationLine(directoryInformation[0]));
    TSeqNo startPos = 0;    // accession separation offset
    TSeqPos offset;         // kmer offset within accession
    TSeqPos kmer_length = primer_cfg.get_kmer_length();
    TKmerLocations::iterator kmer_it = kmer_locations.begin();
    // index to load in order to extract strings
    TIndex index;
    if (!genmap::detail::open(index, seqan::toCString(std::string(io_cfg.get_index_base_path())), seqan::OPEN_RDONLY))
        std::cout << "Error in loading index to index obj.\n", exit(0);
    auto const & text = seqan::indexText(index);
    uint64_t accession_ctr = 0; // fasta header counter
    // zero-based accession counter
    TSeqNo next_accession = seqan::getValueI1<TSeqNo, TSeqPos>(kmer_it->second[0]);
    while (kmer_it != kmer_locations.end())
    {
        auto row = retrieveDirectoryInformationLine(directoryInformation[accession_ctr]);
        // forward id counter and accumulate text offset
        while (accession_ctr != next_accession && accession_ctr < length(directoryInformation))
        {
            ++accession_ctr;
            row = retrieveDirectoryInformationLine(directoryInformation[accession_ctr]);
            startPos += std::get<1>(row); // add sequence length
        }
        if (std::get<0>(row) != fastaFile)
            break;
        assert(accession_ctr == next_accession);

        auto chromosomeNames = std::get<2>(row);
        std::cout << "accession_ctr = " << accession_ctr << ", corresponding to " << chromosomeNames << std::endl;
        std::cout << "next kmer id = " << next_accession << std::endl;
        // sequence internal offset
        offset = seqan::getValueI2<TSeqNo, TSeqPos>(kmer_it->second[0]);
        std::cout << "Access text.concat at " << startPos + offset << " to " << (startPos + offset+kmer_length) << ", text.concat.length = " << length(text.concat) << std::endl;
        auto const & kmer_str = seqan::infixWithLength(text.concat, startPos + offset, kmer_length);
        std::cout << "loc = (" << seqan::getValueI1<TSeqNo, TSeqPos>(kmer_it->second[0]) << ", " << seqan::getValueI2<TSeqNo, TSeqPos>(kmer_it->second[0]) << ") has kmer sequence = " << kmer_str << std::endl;
        // copy kmer into first position of locations vector
        seqan::String<priset::dna> str;
        std::cout << "append kmer sequence to str" << std::endl;
        seqan::append(str, kmer_str);
        std::cout << "assign to 1st position of kmer_locs: " << std::endl;
        // insert into kmer map
        kmer_map[(*kmer_it).first] = TKmer{(*kmer_it).first, str, chemistry::get_Tm(primer_cfg, str)};
        // forward next kmer iterator and abort if no more kmers to resolve
        std::cout << "increment kmer iterator: " << std::endl;
        ++kmer_it;
        std::cout << "query for end: " << std::endl;
        if (kmer_it == kmer_locations.end())
        {
            std::cout << "kmer_locations end reached\n";
            break;
        }
        std::cout << "get next kmer id: " << std::endl;
        assert(kmer_it->second.size() > 0);
        next_accession = seqan::getValueI1<TSeqNo, TSeqPos>(kmer_it->second[0]);
    }
}

// set directory information as needed by genmap's fasta file parsing
void set_directoryInformation(std::string & index_path_base_ids, TDirectoryInformation & directoryInformation)
{
    seqan::open(directoryInformation, seqan::toCString(index_path_base_ids), seqan::OPEN_RDONLY);
    // dummy entry enforces that the mappability is computed for the last file in the while loop.
    seqan::appendValue(directoryInformation, "dummy.entry;0;chromosomename");
}

// split csv row into tokens
void split(std::string const & line, std::vector<std::string> & tokens, std::string const & delimiter = ",")
{
    tokens.resize(0);
    size_t pos = 0, pos_prev = 0;
    std::string token;
    do
    {
        pos = line.find(delimiter, pos_prev);
        if (pos == std::string::npos)
            pos = line.length();
        token = line.substr(pos_prev, pos - pos_prev);
        if (!token.empty())
            tokens.push_back(token);
        pos_prev = pos + delimiter.length();
    }
     while (pos < line.length() && pos_prev < line.length());
}

// create_table helper to build accession ID to accession number map
void create_accID2acc_map(std::unordered_map<TAccID, std::string> & accID2acc, std::unordered_map<TAcc, TAccID> & acc2accID, io_cfg_type const & io_cfg)
{
    std::cout << "create_accID2acc_map\n";
    std::ifstream id_file(io_cfg.get_id_file());
    std::vector<std::string> tokens;
    while (id_file)
    {
        std::string line;
        if (!getline(id_file, line)) break;
        if (line[0] != '#')
        {
            split(line, tokens);
            if (tokens.size() != 2)
                std::cout << "ERROR: unknown id,acc format in " << io_cfg.get_id_file() << std::endl, exit(0);
            TAccID accID = std::stoi(tokens[0]);
            accID2acc[accID] = tokens[1];
            acc2accID[tokens[1]] = accID;
        }
    }
    std::cout << "... done\n";
}

// create_table helper to build accession ID to taxon ID map
void create_accID2taxID_map(std::unordered_map<TAccID, TTaxid> & accID2taxID, std::unordered_set<TTaxid> & taxid_set, std::unordered_map<TAcc, TAccID> const & acc2accID, io_cfg_type const & io_cfg)
{
    std::cout << "create_accID2taxID\n";
    std::ifstream acc_file(io_cfg.get_acc_file());
    std::vector<std::string> tokens;
    // taxid: (ctr_match, ctr_total), ctrs for accessions
    while (acc_file)
    {
        std::string line;
        if (!getline(acc_file, line)) break;
        if (line[0] != '#')
        {
            split(line, tokens);
            TTaxid taxid = std::stoi(tokens[0]);
            taxid_set.insert(taxid);
            for (uint16_t token_idx = 1; token_idx < tokens.size(); ++token_idx)
            {
                TAcc acc = tokens[token_idx];
                accID2taxID[acc2accID.at(acc)] = taxid;
            }
        }
    }
    std::cout << "... done\n";
}

// load taxonomy from file and store as map {taxid: p_taxid}
void create_tax_map(std::unordered_map<TTaxid, TTaxid> & tax_map, io_cfg_type const & io_cfg)
{
    std::cout << "create_tax_map\n";
    std::ifstream tax_file(io_cfg.get_tax_file());
    size_t pos;
    while (tax_file)
    {
        std::string line;
        if (!getline(tax_file, line)) break;
        if (line[0] != '#')
        {
            pos = line.find(",");
            if (pos == std::string::npos)
                continue;
            tax_map[std::stoi(line.substr(0, pos))] = std::stoi(line.substr(pos + 1, std::string::npos));
        }
    }
    std::cout << "... done\n";
}

// write results in csv format
// columns: taxid, fwd, rev, num_IDs_match, num_IDs_total, ID_list
void create_table(io_cfg_type const & io_cfg, TKmerLocations const & kmer_locations, TKmerMap const & kmer_map, TKmerPairs const & pairs)
{
    // TODO: check if unordered_map instead of map
    // load id file for mapping reference IDs (1-based) to accession numbers and vice versa
    std::unordered_map<TAccID, TAcc> accID2acc;
    std::unordered_map<TAcc, TAccID> acc2accID;
    create_accID2acc_map(accID2acc, acc2accID, io_cfg);

    // build dictionary for taxids and counter for assigned accessions
    std::unordered_map<TAccID, TTaxid> accID2taxID;
    std::unordered_set<TTaxid> taxid_set; // taxids with accessions
    create_accID2taxID_map(accID2taxID, taxid_set, acc2accID, io_cfg);

    // load taxonomy as map {taxid: p_taxid}, taxid is root if not in key set
    std::unordered_map<TTaxid, TTaxid> tax_map;
    create_tax_map(tax_map, io_cfg);

    std::cout << "ct1\n";
    // collect single kmer matches for bottom nodes
    std::unordered_map<TKmerID, std::vector<TSeqNo> > kmer2loc; // relates kmer IDs and location IDs
    for (auto it = kmer_locations.begin(); it != kmer_locations.end(); ++it)
    {
        TKmerID kmer_ID = it->first;
        std::vector<TSeqNo> seq_IDs;
        for (TLocation loc : it->second)
            seq_IDs.push_back(seqan::getValueI1<TSeqNo, TSeqPos>(loc));
        kmer2loc[kmer_ID] = seq_IDs;
    }
    std::cout << "ct2\n";
    // 0-based height, correct level info iff a taxonomic node is in the predecessor lineage of another one
    std::unordered_map<TTaxid, uint16_t> leaves; //<taxid, height_from_bottom>

    for (auto const & taxid : taxid_set)
        leaves[taxid] = 0;
    std::cout << "ct3\n";
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

    std::cout << "ct4\n";
    // taxids sorted by level to have correct upstreams stats
    std::vector<std::pair<TTaxid, uint16_t>> leaves_srt_by_level(leaves.size());
    std::copy(leaves.begin(), leaves.end(), leaves_srt_by_level.begin());
    std::sort(leaves_srt_by_level.begin(), leaves_srt_by_level.end(), [](auto const & l1, auto const & l2){ return l1.second < l2.second; });

    using TUpstreamValue = std::pair<uint16_t, uint16_t>; // match_ctr, covered_taxids
    std::unordered_map<TUpstreamKey::THash, TUpstreamValue > upstream_map;

    // temp structure for result collection
    std::vector<TResult> results;

    // get output file
    std::ofstream table;
    table.open(io_cfg.get_table_file());
        std::cout << "ct5\n";
    for (auto const & [taxid, level] : leaves_srt_by_level)
    {
        // write out single primer results
        for (TKmerLocation kmer_location : kmer_locations) // taxid, kmer fixed
        {
            TKmerID kmerID = kmer_location.first;
            // if one of the reference IDs of the kmer's locations is also in tax2accID, then set match counter to 1, else 0
            std::vector<TAcc> acc_by_tax;
            for (TLocation location : kmer_location.second) // loc = seqan::Pair<TSeqNo, TSeqPos>
            {
                TSeqNo accID = seqan::getValueI1<TSeqNo, TSeqPos>(location);
                if (accID2taxID[accID] == taxid)
                    acc_by_tax.push_back(accID2acc[accID]);
            }
            uint16_t match_ctr = acc_by_tax.size() > 0;
            TResult result{taxid, kmer_location.first, 0, match_ctr, 1, acc_by_tax};
            // we may write back accumulated stats if current node is in lineage of already processed, lower-level node
            if (level)
            {
                TUpstreamKey::THash key{TUpstreamKey{taxid, kmer_location.first, 0}.to_string()};
                auto up_it = upstream_map.find(key);
                if (up_it != upstream_map.end())
                {
                    result.match_ctr += up_it->second.first;
                    result.covered_taxids += up_it->second.second;
                }
            }
            results.push_back(result);

            // accumulate stats for upstream until root
            TTaxid taxid_aux = taxid;
//            auto p_it{tax_map.find(taxid_aux)};
            std::cout << "ct6\n";
            while (tax_map.find(taxid_aux) != tax_map.end())
            {
                // proceed with taxonomic parent
                taxid_aux = tax_map[taxid_aux];
                // key for parental stats
                TUpstreamKey::THash key{TUpstreamKey(taxid_aux, kmer_location.first, 0).to_string()};
                auto st_it{upstream_map.find(key)};
                if (st_it != upstream_map.end())
                {

                    st_it->second.first += result.match_ctr;
                    st_it->second.second += result.covered_taxids;
                }
                else
                    upstream_map[key] = TUpstreamValue{result.match_ctr, result.covered_taxids};
            }
            std::cout << "ct7\n";

            // else taxid is root, no further bottom-up accumulation
        }
    }
    std::cout << "ct8\n";
    // all taxids processed, write out leave nodes and upstream stats
    for (TResult result : results){
        std::cout << result.to_string() << std::endl;
        table << result.to_string();
    }
    for (auto const & [key, value] : upstream_map)
    {
        std::cout << key << "," << value.first << "," << value.second << "\n";
        table << key << "," << value.first << "," << value.second << "\n";
    }

    // clear upstream map for new kmer combinations
    upstream_map.clear();
    // collect kmer pair matches for bottom nodes

    // accumulate counters until root is reached
    table.close();

}

}  // namespace priset
