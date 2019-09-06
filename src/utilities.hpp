#pragma once

#include <algorithm>
#include <array>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <set>
#include <stdexcept>
#include <stdint.h>
#include <string>
#include <unordered_set>
#include <vector>

#include <seqan/basic.h>

#include "../submodules/genmap/src/genmap_helper.hpp"

#include "chemistry.hpp"
#include "io_cfg_type.hpp"
#include "types.hpp"

namespace priset
{

// Identify highest set bit without loop
static inline uint64_t log2_asm(uint64_t const x) {
  uint64_t y;
  asm( "\tbsr %1, %0\n"
        : "=r"(y)
        : "r" (x));
  return y;
}

// Execute in terminal and collect command return value.
std::string exec(char const * cmd) {
    std::cout << "Enter util.exec with cmd = " << std::string(cmd) << std::endl;
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe)
        throw std::runtime_error("ERROR: popen() failed!");
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr)
        result += buffer.data();
    std::cout << "util.exec, result = " << result << std::endl;
    return result;
}

// forward declaration
struct primer_cfg_type;
//struct io_cfg_type;

// helper to compute number of unique kmers, insert constant, but not sorted
void unique_kmers(TKmerIDs const & kmerIDs, std::unordered_set<TKmerID> & kmer_set)
{
    kmer_set.clear();
    for (auto kmerID_list : kmerIDs)
        for (TKmerID kmerID : kmerID_list)
            kmer_set.insert(kmerID);
}

// helper to compute number of unique kmers, insert log(n), but sorted
void unique_kmers(TKmerIDs const & kmerIDs, std::set<TKmerID> & kmer_set)
{
    kmer_set.clear();
    for (auto kmerID_list : kmerIDs)
        for (TKmerID kmerID : kmerID_list)
            kmer_set.insert(kmerID);
}

/* Encode a single sequence as a 64 bit integer.
 * Details: encoding schme is \sum_i 4^i*x_i, starting with the first character
 * (little endian) and x_i being the 2 bit representation of 'A' (=0), 'C' (=1),
 * 'G' (=2), and 'T' (=4) ..., 3 = 'G'. E.g., ACGT is encoded as 0*4^0 + 1*4^1 + 2*4^2 + 3*4^2.
 * Non-zero encoded character ('C') is added because of flexible sequence lengths
 * and therefore the necessity to differentiate 'XA' from 'XAA'.
 * Since multiple kmers may start at one position in the reference (which means that they all
 * share the same prefix), the code stores in its 12 highest bits flags for which length
 * the code represents. E.g. if the 1st bit is set, the code represents a string
 * with the length of the shortest possible primer length, if the 5th bit is set, the code also represents
 * a string of minimal primer length plus 4, and so forth.
 *
 * Bit layout:
 * bits [0 .. |seq|-1]  complete sequence
 * bit [|seq|]          closure symbol 'C'
 * bits [60:64]         lower sequence length bound in case of variable length
 */
 uint64_t dna_encoder(TSeq const & seq)
 {
     std::cout << "Enter dna_encoder\n";
     std::cout << "input seq in dna_encoder: " << seq << std::endl;
     uint64_t code = 1ULL << uint64_t(seqan::length(seq) << 1ULL); // stop symbol 'C' = 1
     for (uint64_t i = 0; i < seqan::length(seq); ++i)
     {
         switch (char(seqan::getValue(seq, i))) //char(seq[i]))
         {
             case 'C': code |=  1ULL << (i << 1ULL); break;
             case 'G': code |=  2ULL << (i << 1ULL); break;
             case 'T': code |=  3ULL << (i << 1ULL);
         }
     }

         std::cout << "Return from dna_encoder\n";
     return code;
 }

// return full length sequence, ignore variable length info in leading bits.
TSeq dna_decoder(uint64_t code)
{
    // note that assert converted to nop due to seqan's #define NDEBUG
    if (code == 0ULL)
        throw std::invalid_argument("ERROR: invalid argument for decoder, code > 0.");
    code &= (1ULL << 60ULL) - 1ULL;
    std::array<std::string, 4> sigmas = {"A", "C", "G", "T"};
    TSeq d = "";
    while (code != 1)
    {
        d += sigmas[3 & code];
        code >>= 2;
    }
    return d;
}

void dna_decoder(primer_cfg_type const & primer_cfg, uint64_t code, std::vector<TSeq> & decodes)
{
    // note that assert converted to nop due to seqan's #define NDEBUG
    if (code == 0ULL)
        throw std::invalid_argument("ERROR: invalid argument for decoder, code > 0.");
    decodes.clear();
    uint64_t kmer_length_mask = (code & ~((1ULL << 52ULL) - 1ULL)) >> 52ULL;
    if (code == 0ULL)
        throw std::invalid_argument("ERROR: invalid argument for decoder, code > 0.");
    decodes.clear();
    code &= (1ULL << 52ULL) - 1ULL; // clear leading kmer length information
    TSeq decode = dna_decoder(code); // largest kmer
    uint64_t K = seqan::length(decode); //
    kmer_length_mask >>= primer_cfg.primer_min_length - K;
    while (kmer_length_mask != 0)
    {
        if (kmer_length_mask & 1)
            decodes.push_back(seqan::prefix(decode, K)); //.substr(0, K));
        kmer_length_mask >>= 1;
        --K;
    }
}

uint64_t location_encode(TSeqNo seqNo, TSeqPos seqPos)
{
    return (seqNo << 10ULL) + seqPos;
}


TKmerLength get_kmer_length(uint64_t code)
{
    // note that assert converted to nop due to seqan's #define NDEBUG
    if (code == 0ULL)
        throw std::invalid_argument("ERROR: invalid argument for decoder, code > 0.");
    TKmerLength K = 0;
    assert(code >= 1);
    while (code != 1)
    {
        code >>= 2;
        ++K;
    }
    std::cout << "leave with K = " << K << std::endl;
    return K;
}

// print binary format
template<typename uint_type>
void print_bits(uint_type i)
{
    if (!i)
        std::cout << "0";
    while(i)
        (i & 1) ? std::cout << "1" : std::cout << "0";
    std::cout << std::endl;
}

void print_combinations(primer_cfg_type const & primer_cfg, TKmerIDs const & kmerIDs, TPairs const & pairs) noexcept
{
    std::cout << "KmerID_fwd\t| KmerID_rev\t| Substring Combinations \n";
    std::cout << "----------------------------------------------------\n";
    for (TPair pair : pairs)
    {
        TKmerID kmerID_fwd = std::get<0>(pair);
        TKmerID kmerID_rev = std::get<1>(pair);

        std::cout << TKmerID_fwd << "\t | " << TKmerID_rev << "\t| ";
        TCombinePattern cp = std::get<2>(pair);
        std::vector<std::pair<TKmerLength, TKmerLength>> combinations;
        cp.get_combinations(combinations, primer_cfg.primer_min_length);
        TSeq kmer_fwd = decode(kmerID_fwd);
        TSeq kmer_rev = decode(kmerID_rev);
        for (auto lc : combinations)
            std::cout << "\t\t | \t\t | (kmerID_fwd[" lc.first << ":], kmerID_rev[" << lc.second << "]) = (" << seqan::infixWithLength(kmer_fwd, 0, lc.first) << ", " << seqan::infixWithLength(kmer_rev, 0, lc.second) << ")\n";
    }
}

void print_locations(TLocations & locations)
{
    using key_type = typename TLocations::key_type;
    using TSeqNo = typename seqan::Value<key_type, 1>::Type;
    using TSeqPos = typename seqan::Value<key_type, 2>::Type;
    //using value1_type = typename seqan::Value<typename TLocations::value_type, 1>::Type;
    std::cout << "Locations (SeqNo, SeqPos): [(SeqNo, SeqPos)], [(SeqNo, SeqPos)]:\n";
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

/*
void print_kmer_locations(TKmerLocations const & kmer_locations)
{
    for (TKmerLocations::size_type i = 0; i < kmer_locations.size(); ++i)
    {
        auto kmer_ID = kmer_locations[i].get_kmer_ID();
        std::cout << dna_decoder(kmer_ID) << ": [";
        for (TKmerLocation::size_type j = 0; j < kmer_locations[i].container_size(); ++j)
            std::cout << "(" << kmer_locations[i].accession_ID_at(j) << ", " << kmer_locations[i].kmer_pos_at(j) << ") ";
        std::cout << "]" << std::endl;
    }
}*/

void print_pairs(primer_cfg_type const & primer_cfg, TKmerPairs const & kmer_pairs)
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
        for (TKmerPairs::size_type i = 0; i < (*it).container_size(); ++i)
            std::cout << (*it).accession_ID_at(i) << "\t\t| (" << (*it).kmer_pos_at(i, 1) << ", " << (*it).kmer_pos_at(i, 2) << ")\n";
    }
    if (kmer_pairs.size())
    {
        std::cout << "\n\nkmer ID\t | sequence(s)\n------------------------------\n";
        std::vector<TSeq> decodes;
        for (TKmerID kmer_ID : legend)
        {
            dna_decoder(primer_cfg, kmer_ID, decodes);
            std::cout << kmer_ID << "\t| ";
            for (auto decode : decodes)
                std::cout << decode << " " << std::endl;
        }

    }
}

/*
// Retrieve DNA sequences from txt.concat given a set of locations. Kmer IDs are retrieved from
void lookup_sequences(TKmerLocations & kmer_locations, io_cfg_type const & io_cfg, primer_cfg_type const & primer_cfg)
{
    // load corpus
    seqan::StringSet<seqan::DnaString, seqan::Owner<seqan::ConcatDirect<>>> text;
    fs::path text_path = io_cfg.get_index_txt_path();
    std::cout << "text_path: " << text_path << std::endl;
    seqan::open(text, text_path.string().c_str(), seqan::OPEN_RDONLY);

    TSeqNo kmer_ID;
    TSeqPos kmer_pos;
    TKmerLength K;
    // TODO: remove this later, uniqueness should be already ensured by location construction
    std::unordered_set<TKmerID> seen;
    for (TKmerLocations::iterator kmer_it = kmer_locations.begin(); kmer_it != kmer_locations.end(); ++kmer_it)
    {
        acc_ID = kmer_it->accession_ID_at(0);
        kmer_pos = kmer_it->kmer_pos_at(0);
        K = kmer_it->get_K();
        seqan::DnaString seq = seqan::valueById(text, acc_ID);
        auto const & kmer_str = seqan::infixWithLength(seq, kmer_pos, K);
        uint64_t kmer_code = dna_encoder(kmer_str);
        assert(!seen.contains(kmer_code));
        seen.insert(kmer_code)
        // replace kmer_ID
        kmer_it->set_kmer_ID(kmer_code);
    }
}*/

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
template<typename io_cfg_type>
void create_accID2acc_map(std::unordered_map<TAccID, std::string> & accID2acc, std::unordered_map<TAcc, TAccID> & acc2accID, io_cfg_type const & io_cfg)
{
    //std::cout << "create_accID2acc_map\n";
    std::ifstream id_file(io_cfg.get_id_file());
    std::vector<std::string> tokens;
    std::string line;
    // skip header
    getline(id_file, line);
    while (id_file)
    {
        if (!getline(id_file, line)) break;
        split(line, tokens);
        if (tokens.size() != 2)
            std::cout << "ERROR: unknown id,acc format in " << io_cfg.get_id_file() << std::endl, exit(0);
        TAccID accID = std::stoi(tokens[0]);
        for (size_t i = 1; i < tokens.size(); ++i)
        {
            accID2acc[accID] = tokens[i];
            acc2accID[tokens[i]] = accID;
            //std::cout << "filled both dictionaries with " << accID << " <-> " << tokens[i] << std::endl;
        }
    }
    std::cout << "... done\n";
}

// create_table helper to build accession ID to taxon ID map
template<typename io_cfg_type>
void create_accID2taxID_map(std::unordered_map<TAccID, TTaxid> & accID2taxID, std::unordered_set<TTaxid> & taxid_set, std::unordered_map<TAcc, TAccID> const & acc2accID, io_cfg_type const & io_cfg)
{
    std::cout << "create_accID2taxID\n";
    std::ifstream acc_file(io_cfg.get_acc_file());
    std::cout << "load acc_file: " << io_cfg.get_acc_file() << std::endl;
    std::vector<std::string> tokens;
    // taxid: (ctr_match, ctr_total), ctrs for accessions
    std::string line;
    // skip header
    getline(acc_file, line);
    while (acc_file)
    {
        if (!getline(acc_file, line)) break;

        split(line, tokens);
        TTaxid taxid = std::stoi(tokens[0]);
        //std::cout << "taxid = " << taxid << std::endl;
        taxid_set.insert(taxid);
        for (uint16_t token_idx = 1; token_idx < tokens.size(); ++token_idx)
        {
            TAcc acc = tokens[token_idx];
            //std::cout << "acc = " << acc << std::endl;
            // TODO: observation - there are accessions (without version suffix) that do not have fasta entries in DB
            if (acc2accID.find(acc) == acc2accID.end())
                continue; //std::cout << "ERROR: accession " << acc << " not in acc2accID dictionary!" << std::endl, exit(0);
            accID2taxID[acc2accID.at(acc)] = taxid;
        }
    }
    std::cout << "... done\n";
}

// load taxonomy from file and store as map {taxid: p_taxid}
template<typename io_cfg_type>
void create_tax_map(std::unordered_map<TTaxid, TTaxid> & tax_map, io_cfg_type const & io_cfg)
{
    std::cout << "create_tax_map\n";
    std::ifstream tax_file(io_cfg.get_tax_file());
    size_t pos;
    std::string line;
    // skip header line
    getline(tax_file, line);
    while (tax_file)
    {
        if (!getline(tax_file, line)) break;
        pos = line.find(",");
        if (pos == std::string::npos)
            continue;
        tax_map[std::stoi(line.substr(0, pos))] = std::stoi(line.substr(pos + 1, std::string::npos));
    }
    std::cout << "... done\n";
}

// accumulate statistics upstream for both container types - TKmerLocations and TKmerPairs
template<typename TKmerContainer, typename io_cfg_type>
void accumulation_loop(TKmerContainer const & kmer_container, std::vector<std::pair<TTaxid, uint16_t>> const & leaves, std::unordered_map<TTaxid, TTaxid> const & tax_map, std::unordered_map<TAccID, TTaxid> const & accID2taxID, std::unordered_map<TAccID, TAcc> const & accID2acc, io_cfg_type const & io_cfg)
{
    if (!kmer_container.size())
        return;
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
        for (typename TKmerContainer::value_type kmer_location : kmer_container) // taxid, kmer fixed
        {
            // CHECK: behaviour correct for kmer_locations with only one kmer_ID
            TKmerID kmerID1 = kmer_location.get_kmer_ID1();
            TKmerID kmerID2 = kmer_location.get_kmer_ID2();

            // collect all accessions (not accession IDs) assigned to current taxid where kmer matches
            std::vector<TAcc> acc_by_tax;
            for (typename TKmerLocation::size_type i = 0; i < kmer_location.container_size(); ++i)  //TLocation location : kmer_location.second) // loc = seqan::Pair<TSeqNo, TSeqPos>
            {
                TSeqNo accID = kmer_location.accession_ID_at(i);  //seqan::getValueI1<TSeqNo, TSeqPos>(location);
                if (accID2taxID.at(accID) == taxid)
                    acc_by_tax.push_back(accID2acc.at(accID));
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
            while (tax_map.find(taxid_aux) != tax_map.end())
            {
                // proceed with taxonomic parent
                taxid_aux = tax_map.at(taxid_aux);
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
    table.open(io_cfg.get_result_file(), std::ios_base::app);

    // flush leave node results
    for (TResult result : results){
        //std::cout << result.to_string() << std::endl;
        table << result.to_string();
    }
    // flush inner node results
    for (auto const & [key, value] : upstream_map)
    {
        //std::cout << key << "," << value.first << "," << value.second << "\n";
        table << key << "," << value.first << "," << value.second << "\n";
    }

    table.close();
}


// Result output helper for writing primer infos.
template<typename io_cfg_type, typename primer_cfg_type>
void write_primer_info_file(io_cfg_type const & io_cfg, primer_cfg_type const & primer_cfg, TKmerIDs const & kmerIDs)
{
    std::ofstream primer_table;
    primer_table.open(io_cfg.get_primer_info_file());
    primer_table << "kmer_ID,kmer_sequence,Tm\n";
    std::set<TKmerID> kmer_ordered_set;
    unique_kmers(kmerIDs, kmer_ordered_set);
    uint64_t i = 0;
    for (TKmerID kmerID : kmer_ordered_set)
    {
        primer_table << (i++) << "," << dna_decoder(kmerID) << "," << get_Tm(primer_cfg, kmerID) << "\n";
    }
    primer_table.close();
    std::cout << "STATUS: primer_info.csv written to\t" << io_cfg.get_primer_info_file() << std::endl;
}

/*
    Write result table with columns: taxid, fwd, rev, matches, coverage, ID_list and
    primer info file with columns kmer_id (1-based), sequence and melting temperature.
*/
template<typename io_cfg_type, typename primer_cfg_type>
void create_table(io_cfg_type const & io_cfg, primer_cfg_type const & primer_cfg, TReferences const & references, TKmerIDs const & kmerIDs, TPairs const & pairs)
{
    std::set<TKmerID> kmer_ordered_set;
    unique_kmers(kmerIDs, kmer_ordered_set);
    std::cout << "STATUS: will write " << kmer_ordered_set.size() << " single primer results and " << pairs.size() << " primer pair results\n";
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
    accumulation_loop<TKmerPairs>(kmer_pairs, leaves_srt_by_level, tax_map, accID2taxID, accID2acc, io_cfg);
    std::cout << "STATUS: table.csv written to\t" << io_cfg.get_result_file() << std::endl;

    // write primer info file
    write_primer_info_file(io_cfg, primer_cfg, kmer_locations);
    */
}

}  // namespace priset
