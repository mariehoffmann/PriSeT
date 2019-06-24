#pragma once

#include <algorithm>
#include <array>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <seqan/basic.h>

#include "../submodules/genmap/src/genmap_helper.hpp"

#include "chemistry.hpp"
#include "io_cfg_type.hpp"
#include "types.hpp"

namespace priset
{

// Execute in terminal and collect command return value.
std::string exec(char const * cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe)
        throw std::runtime_error("ERROR: popen() failed!");
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr)
        result += buffer.data();
    return result;
}

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
    for (TKmerLocations::size_type i = 0; i < kmer_locations.size(); ++i)
    {
        const TKmer kmer{kmer_map.at(kmer_locations[i].get_kmer_ID())};
        std::cout << kmer.seq << ": [";
        for (TKmerLocation::size_type j = 0; j < kmer_locations[i].container_size(); ++j)
            std::cout << "(" << kmer_locations[i].accession_ID_at(j) << ", " << kmer_locations[i].kmer_pos_at(j) << ") ";

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
        std::cout << "(" << (*it).get_kmer_ID1() << ", " << (*it).get_kmer_ID2() << ")\t\t| ";
        legend.insert((*it).get_kmer_ID1());
        legend.insert((*it).get_kmer_ID2());
        for (auto i = 0; i < (*it).container_size(); ++i)
            std::cout << (*it).accession_ID_at(i) << "\t\t| (" << (*it).kmer_pos_at(i, 1) << ", " << (*it).kmer_pos_at(i, 2) << ")\n";
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
    TSeqNo next_accession = kmer_it->accession_ID_at(0); //seqan::getValueI1<TSeqNo, TSeqPos>(kmer_it->second[0]);
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
        //offset = seqan::getValueI2<TSeqNo, TSeqPos>(kmer_it->second[0]);
        offset = kmer_it->kmer_pos_at(0);
        std::cout << "Access text.concat at " << startPos + offset << " to " << (startPos + offset + kmer_length) << ", text.concat.length = " << length(text.concat) << std::endl;
        auto const & kmer_str = seqan::infixWithLength(text.concat, startPos + offset, kmer_length);
        std::cout << "loc = (" << kmer_it->accession_ID_at(0) << ", " << offset << ") has kmer sequence = " << kmer_str << std::endl;
        // copy kmer into first position of locations vector
        seqan::String<priset::dna> str;
        std::cout << "append kmer sequence to str" << std::endl;
        seqan::append(str, kmer_str);
        std::cout << "assign to 1st position of kmer_locs: " << std::endl;
        // insert into kmer map
//        kmer_map[(*kmer_it).first] = TKmer{(*kmer_it).first, str, chemistry::get_Tm(primer_cfg, str)};
        kmer_map[kmer_it->get_kmer_ID()] = TKmer{kmer_it->get_kmer_ID(), str, chemistry::get_Tm(primer_cfg, str)};

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
        assert(kmer_it->container_size() > 0);
        next_accession = kmer_it->accession_ID_at(0); //seqan::getValueI1<TSeqNo, TSeqPos>(kmer_it->second[0]);
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

// accumulate statistics upstream for both container types - TKmerLocations and TKmerPairs
template<typename TKmerContainer>
void accumulation_loop(TKmerContainer const & kmer_container, std::vector<std::pair<TTaxid, uint16_t>> const & leaves, std::unordered_map<TTaxid, TTaxid> & tax_map, io_cfg_type const & io_cfg)
{
    // type for upstream stats collection: (match_ctr, covered_taxids)
    using TUpstreamValue = std::pair<uint16_t, uint16_t>;
    // for each taxid, primer_fwd, primer_rev (as string) combination store counters for matches and coverage
    std::unordered_map<TUpstreamKey::THash, TUpstreamValue > upstream_map;
    // temp structure for result row collection
    std::vector<TResult> results;

    for (auto const & [taxid, level] : leaves)
    {
        // write out single primer results
        // value_type is either TKmerLocation or TKmerPair
        for (TKmerContainer::value_type kmer_location : kmer_container) // taxid, kmer fixed
        {

            TKmerID kmerID1 = kmer_location.get_kmer_ID1;
            TKmerID kmerID2 = kmer_location.get_kmer_ID2;

            // collect all accessions (not accession IDs) assigned to current taxid where kmer matches
            std::vector<TAcc> acc_by_tax;
            for (kmer_location::size_type i = 0; i < kmer_location.container_size(); ++i)  //TLocation location : kmer_location.second) // loc = seqan::Pair<TSeqNo, TSeqPos>
            {
                TSeqNo accID = kmer_location.accession_ID_at(i);  //seqan::getValueI1<TSeqNo, TSeqPos>(location);
                if (accID2taxID[accID] == taxid)
                    acc_by_tax.push_back(accID2acc[accID]);
            }
            uint16_t match_ctr = acc_by_tax.size() > 0;
            TResult result{taxid, kmerID1, kmerID2, match_ctr, 1, acc_by_tax};
            // we may write back accumulated stats if current node is in lineage of already processed, lower-level node
            if (level)
            {
                TUpstreamKey::THash key{TUpstreamKey{taxid, kmerID1, kmerID2}.to_string()};
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
                TUpstreamKey::THash key{TUpstreamKey(taxid_aux, kmerID1, kmerID2).to_string()};
                auto st_it{upstream_map.find(key)};
                if (st_it != upstream_map.end())
                {
                    st_it->second.first += result.match_ctr;
                    st_it->second.second += result.covered_taxids;
                }
                else
                    upstream_map[key] = TUpstreamValue{result.match_ctr, result.covered_taxids};
            }
            // else taxid is root, no further bottom-up accumulation
        }
    }

    // create output stream to result table and append
    std::ofstream table;
    table.open(io_cfg.get_table_file(), ios_base::ate);

    // flush leave node results
    for (TResult result : results){
        std::cout << result.to_string() << std::endl;
        table << result.to_string();
    }
    // flush inner node results
    for (auto const & [key, value] : upstream_map)
    {
        std::cout << key << "," << value.first << "," << value.second << "\n";
        table << key << "," << value.first << "," << value.second << "\n";
    }

    table.close();

}

// write results in csv format
// columns: taxid, fwd, rev, matches, coverage, ID_list
// TODO: write out kmer IDs and DNA sequences
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

    std::cout << "ct8\n";
    // write result table header
    std::ofstream table;
    table.open(io_cfg.get_table_file());
    table << "#taxid, fwd, rev, matches, coverage, ID_list\n";
    table.close();

    // collect single kmer matches
    accumulation_loop<TKmerLocations>(kmer_locations, leaves_srt_by_level, tax_map, io_cfg);

    // collect kmer pair matches for bottom nodes
    accumulation_loop<TKmerPairs>(kmer_pairs, leaves_srt_by_level, tax_map, io_cfg);

}

}  // namespace priset
