#pragma once

#include <iostream>

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

void print_kmer_locations(TKmerLocations const & kmer_locations, TKmerMap & kmer_map)
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
void split(std::string const & line, std::string const & delimiter = ",", std::vector<std::string> & tokens)
{
    tokens.resize(0);
    size_t pos = 0, pos_old = 0;
    std::string token;
    do
    {
        pos = line.find(delimiter, pos_old);
        if (pos == string::npos)
            pos = line.length();
        token = line.substr(pos_old, pos - pos_old);
        if (!token.empty())
            tokens.push_back(token);
        pos_prev = pos + delim.length();
    }
     while (pos < line.length() && pos_old < line.length());
}

// Filter accessions from location vector by given taxid and transform TAccID to TAcc.
void filter_acc_by_tax(std::vector<TLocation> const & locations, std::unordered_map<TTaxid, std::vector<TAccID > > const & tax2accID, TTaxid const & taxid, std::unordered_map<TAccID, std::string> const & accID2acc, std::vector<std::string> & result)
{
    result.resize(0);
    // assert non-monotonouos TAccID sequence
    if (!std::is_sorted(tax2accID[taxid]))
        std::cout << "ERROR: accesion ID list in map not sorted!\n", exit(0);
    for (auto it = locations.begin() + 1; it != locations.end(); ++it)
        if (seqan::getValueI1<TSeqNo, TSeqPos>(*(it-1)) > seqan::getValueI1<TSeqNo, TSeqPos>(*it)))
            std::cout << "ERROR: accesion ID list in TLocations not sorted!\n", exit(0);

    auto it1 = locations.begin();
    auto it2 = tax2accID.at(taxid).begin();
    while (it1 != locations.end() && it2 != tax2accID.at(taxid).end())
    {
        if (seqan::getValueI1<TSeqNo, TSeqPos>(*it1) == *it2)
            result.push_back(accID2acc.at(*it2));
        while (it1 != locations.end() && seqan::getValueI1<TSeqNo, TSeqPos>(*it1) <
                *it2) ++it1;
        while (it2 != tax2accID.at(taxid).end() && (*it2) <
                seqan::getValueI1<TSeqNo, TSeqPos>(*it1)) ++it2;
    }
}

// write results in csv format
// columns: taxid, fwd, rev, num_IDs_match, num_IDs_total, ID_list
void create_table(io_cfg_type const & io_cfg, TKmerLocations const & kmer_locations, TKmerMap const & kmer_map, TPairs const & pairs)
{
    // load id file for mapping reference IDs (1-based) to accession numbers
    std::unordered_map<std::string, TAccID> acc2id;
    std::unordered_map<TAccID, std::string> accID2acc;
    ifstream id_file(io_cfg.get_id_file());
    std::vector<std::string> tokens;
    while (id_file)
    {
        std::string line;
        if (!getline(id_file, line)) break;
        if (line[0] != '#')
        {
            split(line, tokens);
            if (tokens.size() != 2)
                std::cout << "ERROR: unknown id,acc format in " << id_file << std::endl, exit(0);
            acc2id[atoi(tokens[0])] = tokens[1];
            acc2id_rev[tokens[1]] = atoi(tokens[0]);
        }
    }

    // build dictionary for taxids and counter for assigned accessions
    ifstream acc_file(io_cfg.get_acc_file());
    std::vector<std::string> tokens;
    // taxid: (ctr_match, ctr_total), ctrs for accessions
    std::unordered_map<TTaxid, std::vector<TAccID > > tax2accID;
    while (acc_file)
    {
        std::string line;
        if (!getline(acc_file, line)) break;
        if (line[0] != '#')
        {
            split(line, tokens);
            auto taxid = atoi(tokens[0]);
            std::vector<TAccID> accessionIDs(tokens.size() - 1);
            tax2accID[taxid] = std::transform(tokens.begin() + 1, tokens.end(), accessionIDs.begin(), [](std::string const & acc){ return acc2id[acc];});
        }
    }

    // load taxonomy as map {taxid: p_taxid}
    ifstream tax_file(io_cfg.get_tax_file());
    std::unordered_map<TTaxid, TTaxid> tax_map;
    size_t pos;
    while (tax_file)
    {
        std::string line;
        if (!getline(tax_file, line)) break;
        pos = line.find(",")
        if (pos == string::npos)
            continue;
        tax_map[atoi(line.substr(0, pos))] = atoi(line.substr(pos + 1, string::npos))
    }

    // get output file
    std::ofstream table;
    table.open(io_cfg.get_table_path());
    // collect single kmer matches for bottom nodes
    std::unordered_map<TKmerID, std::vector<TSeqNo> > kmer2loc; // relates kmer IDs and location IDs
    for (auto it = kmer_locations.begin(); it != kmer_locations.end(); ++it)
    {
        TKmerID kmer_ID = it.first;
        std::vector<TSeqNo> seq_IDs;
        for (auto it_loc : it.second)
            seq_IDs.push_back(it_loc.first);
        kmer2loc[kmer_ID] = seq_IDs;
    }

    // collect
    std::vector<uint16_t> stack, stack_new;
    for (auto const & [tax, accs] : tax2acc)
        stack.push_back(tax);
    std::sort(stack.begin(), stack.end());
    while (stack.size())
    {
        for (auto taxid : stack)
        {
            std::vector<TResult> results;
            // write out single primer results
            for (TKmerLocation kmer_loc : kmer_locations) // taxid, kmer fixed
            {
                TKmerID kmerID = kmer_loc.first;
                // if one of the reference IDs of the kmer's locations is also in tax2accID, then set match counter to 1, else 0
                std::vector<std::string> acc_by_tax;
                filter_acc_by_tax(kmer_loc.second, tax2accID, taxid, accID2acc, joint_accs);

                result.push_back(TResult{taxid, kmer_loc.first, acc_by_tax.size() == 1, 1, acc_by_tax});
            }
            // replace nodes by their parent taxids
            // TODO: continue here
        }

    }
    table << ;




    // accumulate counters until root is reached

    // collect kmer pair matches for bottom nodes

    // accumulate counters until root is reached
    table.close()
}

}  // namespace priset
