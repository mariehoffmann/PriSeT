#pragma once

#include <algorithm>
#include <array>
#include <bitset>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <stdint.h>
#include <string>
#include <unistd.h>
#include <unordered_set>
#include <unordered_map>
#include <vector>

#include <seqan/basic.h>

#include "../submodules/genmap/src/genmap_helper.hpp"

#include "chemistry.hpp"
#include "types/all.hpp"

// Split prefix and code given a kmerID.
#define split_kmerID(kmerID) std::pair<uint64_t, uint64_t>{(uint64_t)kmerID & PREFIX_SELECTOR, (uint64_t)kmerID & ~PREFIX_SELECTOR}

// The longest encoded k-mer length expressed in 2 bit format, i.e. 16 bp are 32 bits.
#define encoded_length(kmerID) WORD_SIZE - 1 - __builtin_clzl(kmerID & ~PREFIX_SELECTOR);

namespace priset
{

// Identify highest set bit without loop
// TODO: cmp number of CPU cycles with fct for leading zeros (+1)
/*
static inline uint64_t log2_asm(uint64_t const x) {
  uint64_t y;
  asm( "\tbsr %1, %0\n"
        : "=r"(y)
        : "r" (x));
  return y;
}
*/

std::string kmerID2str(TKmerID kmerID);

// extern inline std::pair<uint64_t, uint64_t> split(TKmerID);

// helper: delete lengths bits including the one representing 2bit encode l and larger ones
extern inline void delete_length_bits(TKmerID & kmerID, uint8_t l)
{
    // std::cout << "call delete_length_bits with l = " << int(l) << std::endl;
    auto [prefix, code] = split_kmerID(kmerID);
    if (l <= (PRIMER_MIN_LEN << 1))
        kmerID = code;
    else
    {
        uint64_t offset = PRIMER_MAX_LEN + 1 - (l >> 1);
        kmerID = (prefix & (PREFIX_SELECTOR << offset)) | code;
    }
    // std::cout << "delete_length_bits: " << kmerID2str(kmerID) << std::endl;
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

// Trim to true length
extern inline void trim_to_true_length(TKmerID & kmerID)
{
    // std::cout << "trim2tl input: ";
    // std::cout << kmerID2str(kmerID) << std::endl;
    auto [prefix, code] = split_kmerID(kmerID);
    if (!prefix)
        return;
    auto l_max = (26 - ffsll(prefix >> 54)) << 1;
    auto l_enc = WORD_SIZE - __builtin_clzll(code) - 1;
    if (l_max > l_enc || l_enc < int(PRIMER_MIN_LEN << 1))
        throw std::runtime_error("ERROR: Code shorter than 16 or largest encoded length!");
    code >>= (l_enc - l_max);
    kmerID = prefix | code;
}

// extern inline std::pair<uint64_t, uint64_t> split(TKmerID kmerID)
// {
//     //std::cout << "split: " << bits2str(kmerID) << std::endl;
//     return std::pair<uint64_t, uint64_t>{kmerID & PREFIX_SELECTOR, kmerID & ~PREFIX_SELECTOR};
// }

// forward declaration
struct PrimerConfig;

// TODO: move to dna.hpp

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
std::string dna_decoder(uint64_t code, uint64_t const mask);

uint64_t dna_encoder(seqan::String<priset::dna> const & seq)
{
    uint64_t code(0);
     for (uint64_t i = 0; i < seqan::length(seq); ++i)
     {
         switch (char(seqan::getValue(seq, i)))
         {
             case 'C': code |= 1ULL; break;
             case 'G': code |= 2ULL; break;
             case 'T': code |= 3ULL;
         }
         code <<= 2;

     }
     code >>= 2;
     code |= 1ULL << uint64_t(seqan::length(seq) << 1); // stop symbol 'C' = 1
     return code;
}

// with length bit in prefix
uint64_t dna_encoder_with_lbit(seqan::String<priset::dna> const & seq)
{
    uint64_t code(0);
    auto l = seqan::length(seq);
    for (uint64_t i = 0; i < l; ++i)
    {
         switch (char(seqan::getValue(seq, i)))
         {
             case 'C': code |= 1ULL; break;
             case 'G': code |= 2ULL; break;
             case 'T': code |= 3ULL;
         }
         code <<= 2;
     }
     // revert last shift
     code >>= 2;
     // add length bit and stop symbol 'C' = 1
     code |= 1ULL << (63 - (l - PRIMER_MIN_LEN)) | (1ULL << uint64_t(l << 1));
     return code;
}

// Correct encoded length to length indicated by the mask and remove length information.
// If no mask is given (mask = 0), the code is truncated to its largest length encoded in head,
// otherwise it's truncated to the largest length given by the mask (same format like code head).
// A leading bit ('C') remains in the code to signal the end.
extern inline uint64_t get_code(uint64_t const kmerID, uint64_t mask = 0)
{
    auto [prefix, code] = split_kmerID(kmerID);
    if (!code)
    {
        std::cout << "ERROR: expected kmerID not zero\n";
        exit(0);
    }
    uint8_t enc_l = (WORD_SIZE - 1 - __builtin_clzl(code)) >> 1; // encoded length
    if (!mask)
        mask = ONE_LSHIFT_63 >> (enc_l - PRIMER_MIN_LEN);
    uint8_t mask_l = __builtin_clzl(mask) + PRIMER_MIN_LEN;      // selected length
    return (enc_l == mask_l) ? code : code >> ((enc_l - mask_l) << 1);   // kmer length correction
}

// return full length sequence, ignore variable length info in leading bits if present.
// Note: kmer code is only length trimmed when a mask is given.
std::string dna_decoder(uint64_t const code_, uint64_t const mask = 0)
{
    // note that assert converted to nop due to seqan's #define NDEBUG
    if (code_ == 0ULL)
        throw std::invalid_argument("ERROR: invalid argument for decoder, code == 0.");
    uint64_t code = code_;
    if (mask)
        code = get_code(code_, mask);
    else
        code &= ~PREFIX_SELECTOR;

    // TODO: use global
    std::array<char, 4> alphabet = {'A', 'C', 'G', 'T'};
    uint8_t n = (63 - __builtin_clzl(code)) >> 1;
    char seq[n];
    for (uint8_t i = 1; i <= n; ++i, code >>= 2)
        seq[n - i] = alphabet[3 & code];
    return std::string(seq, n);
}

std::string kmerID2str(TKmerID kmerID)
{
    return bits2str(kmerID >> 54) + "|" + dna_decoder(kmerID);
}

extern inline void dna_decoder(uint64_t kmerID, std::vector<TSeq> & decodes)
{
    // note that assert converted to nop due to seqan's #define NDEBUG
    if (kmerID == 0ULL)
        throw std::invalid_argument("ERROR: invalid argument for decoder, code > 0.");
    decodes.clear();
    uint64_t kmer_length_mask = (kmerID & ~((1ULL << 52ULL) - 1ULL)) >> 52ULL;
    if (kmerID == 0ULL)
        throw std::invalid_argument("ERROR: invalid argument for decoder, code > 0.");
    decodes.clear();
    uint64_t code = kmerID & ~PREFIX_SELECTOR; // clear leading kmer length information
    TSeq decode = dna_decoder(code); // largest kmer
    uint64_t K = seqan::length(decode); //
    kmer_length_mask >>= PRIMER_MIN_LEN - K;
    while (kmer_length_mask != 0)
    {
        if (kmer_length_mask & 1)
            decodes.push_back(seqan::prefix(decode, K)); //.substr(0, K));
        kmer_length_mask >>= 1;
        --K;
    }
}

// Assume that sequence IDs and position indices do not exceed 32 bits.
extern inline uint64_t location_encode(TSeqNo seqNo, TSeqPos seqPos)
{
    return (seqNo << 32) | seqPos;
}

// print binary format
template<typename uint_type>
std::string bits2str(uint_type i)
{
    std::string s = "";
    if (!i)
        return "0";
    while(i)
    {
        s += (i & 1) ? "1" : "0";
        i >>= 1;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

template<typename PairList, typename TKmerIDs>
void print_combinations(TKmerIDs const & kmerIDs, PairList const & pairs) noexcept
{
    std::cout << "Reference ID\t| KmerID_fwd\t| KmerID_rev\t| Substring Combinations \n";
    std::cout << "----------------------------------------------------------------------\n";
    std::cout << "pairs.size = " << pairs.size() << std::endl;
    for (auto pair : pairs)
    {
        TKmerID kmerID_fwd = kmerIDs.at(pair.reference).at(pair.r_fwd - 1);
        TKmerID kmerID_rev = kmerIDs.at(pair.reference).at(pair.r_rev - 1);
        std::cout << pair.reference << "\t | " << kmerID_fwd << "\t | " << kmerID_rev << "\t| ";
        std::vector<std::pair<uint8_t, uint8_t>> combinations;
        pair.cp.get_combinations(combinations);
        TSeq kmer_fwd = dna_decoder(kmerID_fwd);
        TSeq kmer_rev = dna_decoder(kmerID_rev);

        for (auto lc : combinations)
        {
            std::cout << "\t | \t\t\t |  \t\t\t | (kmerID_fwd[" << lc.first << ":], kmerID_rev[" << lc.second << "]) = (" << seqan::infixWithLength(kmer_fwd, 0, lc.first + PRIMER_MIN_LEN) << ", " << seqan::infixWithLength(kmer_rev, 0, lc.second + PRIMER_MIN_LEN) << ")\n";
        }
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
void print_kmer_locations(KmerLocations const & kmer_locations)
{
    for (KmerLocations::size_type i = 0; i < kmer_locations.size(); ++i)
    {
        auto kmer_ID = kmer_locations[i].get_kmer_ID();
        std::cout << dna_decoder(kmer_ID) << ": [";
        for (KmerLocation::size_type j = 0; j < kmer_locations[i].container_size(); ++j)
            std::cout << "(" << kmer_locations[i].accession_ID_at(j) << ", " << kmer_locations[i].kmer_pos_at(j) << ") ";
        std::cout << "]" << std::endl;
    }
}*/

/*
// Retrieve DNA sequences from txt.concat given a set of locations. Kmer IDs are retrieved from
void lookup_sequences(KmerLocations & kmer_locations, IOConfig const & io_cfg, PrimerConfig const & primer_cfg)
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
    for (KmerLocations::iterator kmer_it = kmer_locations.begin(); kmer_it != kmer_locations.end(); ++kmer_it)
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
template<typename IOConfig>
void create_accID2acc_map(std::unordered_map<AccessionID, std::string> & accID2acc, std::unordered_map<Accession, AccessionID> & acc2accID, IOConfig const & io_cfg)
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
        AccessionID accID = std::stoi(tokens[0]);
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
template<typename IOConfig>
void create_accID2taxID_map(std::unordered_map<AccessionID, Taxid> & accID2taxID, std::unordered_set<Taxid> & taxid_set, std::unordered_map<Accession, AccessionID> const & acc2accID, IOConfig const & io_cfg)
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
        Taxid taxid = std::stoi(tokens[0]);
        //std::cout << "taxid = " << taxid << std::endl;
        taxid_set.insert(taxid);
        for (uint16_t token_idx = 1; token_idx < tokens.size(); ++token_idx)
        {
            Accession acc = tokens[token_idx];
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
template<typename IOConfig>
void create_tax_map(std::unordered_map<Taxid, Taxid> & tax_map, IOConfig const & io_cfg)
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

// accumulate statistics upstream for both container types - KmerLocations and TKmerPairs
template<typename TKmerContainer, typename IOConfig, typename Result>
void accumulation_loop(TKmerContainer const & kmer_container, std::vector<std::pair<Taxid, uint16_t>> const & leaves, std::unordered_map<Taxid, Taxid> const & tax_map, std::unordered_map<AccessionID, Taxid> const & accID2taxID, std::unordered_map<AccessionID, Accession> const & accID2acc, IOConfig const & io_cfg)
{
    if (!kmer_container.size())
        return;
    // type for upstream stats collection: (match_ctr, covered_taxids)
    using TUpstreamValue = std::pair<uint16_t, uint16_t>;
    // for each taxid, primer_fwd, primer_rev (as string) combination store counters for matches and coverage
    std::unordered_map<TUpstreamKey::THash, TUpstreamValue > upstream_map;
    // temp structure for result row collection
    std::vector<Result> results;

    for (auto const & [taxid, level] : leaves)
    {
        // write out single primer results
        // value_type is either KmerLocation or TKmerPair
        for (typename TKmerContainer::value_type kmer_location : kmer_container) // taxid, kmer fixed
        {
            // CHECK: behaviour correct for kmer_locations with only one kmer_ID
            TKmerID kmerID1 = kmer_location.get_kmer_ID1();
            TKmerID kmerID2 = kmer_location.get_kmer_ID2();

            // collect all accessions (not accession IDs) assigned to current taxid where kmer matches
            std::vector<Accession> acc_by_tax;
            for (typename KmerLocation::size_type i = 0; i < kmer_location.container_size(); ++i)  //TLocation location : kmer_location.second) // loc = seqan::Pair<TSeqNo, TSeqPos>
            {
                TSeqNo accID = kmer_location.accession_ID_at(i);  //seqan::getValueI1<TSeqNo, TSeqPos>(location);
                if (accID2taxID.at(accID) == taxid)
                    acc_by_tax.push_back(accID2acc.at(accID));
            }
            uint16_t match_ctr = acc_by_tax.size() > 0;
            Result result{taxid, kmerID1, kmerID2, match_ctr, 1, acc_by_tax};
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
            Taxid taxid_aux = taxid;
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
    for (Result result : results){
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

// Build sorted set of unique encoded kmer IDs (unpacked and with null header),
// such that real length is indicated by highest set '1' add an even bit position
// between [2*primer_min_length : 2 : 2*primer_max_length].
template<typename TKmerIDs>
void unique_kmers(TKmerIDs const & kmerIDs, std::set<uint64_t> & set)
{
    //std::vector<std::deque<TKmerID>> TKmerIDs;
    for (auto it_ref = kmerIDs.begin(); it_ref < kmerIDs.end(); ++it_ref)
    {
        for (auto it_kmerID = it_ref->begin(); it_kmerID <= it_ref->end(); ++it_kmerID)
        {
            auto [prefix, code] = split_kmerID(*it_kmerID);
            bool start_shift = false; // since kmer IDs represent only the longest kmer they encode and not necessarily the longest possible primer length, we start truncating the ID after we have seen the first length bit
            while (prefix)
            {
                if (prefix & 1)
                {
                    start_shift = true;
                    set.insert(code);
                }
                prefix >>= 1;
                if (start_shift)
                    code >>= 2;
            }
        }
    }
}

// Count non-unique kmers collected for all references.
uint64_t get_num_kmers(TKmerIDs const & kmerIDs)
{
    uint64_t ctr = 0;
    for (auto kmerIDs_per_reference : kmerIDs)
    {
        for (auto kmerID : kmerIDs_per_reference)
            ctr += __builtin_popcountll(kmerID >> 54);
    }
    return ctr;
}

// Count pairs, i.e. sum of all combination pattern bits of each pair.
template<typename PairList>
uint64_t get_num_pairs(PairList const & pairs)
{
    uint64_t ctr = 0;
    for (auto pair : pairs)
        ctr += pair.cp.size();
    return ctr;
}

// i,j are codes with no prefixes
extern inline uint64_t hash_pair(uint64_t i, uint64_t j)
{
    return (i << __builtin_clzll(i + 1)) ^ j;
}

// Collect hash values of unique pairs and their frequencies.
template<typename PairList, typename TKmerIDs, typename TKmerLength>
void unique_pairs(PairList const & pairs, TKmerIDs const & kmerIDs, std::unordered_map<uint64_t, uint32_t> & code_pairs)
{
    for (auto pair : pairs)
    {
        std::vector<std::pair<uint8_t, uint8_t>> combinations;
        pair.cp.get_combinations(combinations);
        uint64_t code_fwd, code_rev;
        for (std::pair<TKmerLength, TKmerLength> c : combinations)
        {
        //    std::cout << "pair.reference = " << pair.reference << ", size kmerIDs of this reference = " << kmerIDs.at(pair.reference).size();
        //    std::cout << ", try to access r_fwd = " << pair.r_fwd - 1 << ", r_rev = " << pair.r_rev - 1 << std::endl;
            code_fwd = get_code(kmerIDs.at(pair.reference).at(pair.r_fwd - 1), ONE_LSHIFT_63 >> c.first);
            code_rev = get_code(kmerIDs.at(pair.reference).at(pair.r_rev - 1), ONE_LSHIFT_63 >> c.second);
            uint64_t key = hash_pair(code_fwd, code_rev);
            if (code_pairs.find(key) != code_pairs.end())
                ++code_pairs[key];
            else
                code_pairs[key] = 1;
        }
    }
}

// Accumulate unique pair frequencies.
template<typename PairList, typename TKmerIDs, typename TKmerLength>
void count_unique_pairs(PairList const & pairs, TKmerIDs const & kmerIDs)
{
    std::unordered_map<uint64_t, uint32_t> code_pairs;
    unique_pairs<PairList, TKmerIDs, TKmerLength>(pairs, kmerIDs, code_pairs);
    // count frequencies
    std::map<uint64_t, uint64_t> freq_ctrs;
    for (auto it = code_pairs.begin(); it != code_pairs.end(); ++it)
    {
        if (freq_ctrs.find(it->second) != freq_ctrs.end())
            ++freq_ctrs[it->second];
        else
            freq_ctrs[it->second] = 1;
    }
    std::cout << "frequency:\t";
    for (auto it = freq_ctrs.begin(); it != freq_ctrs.end(); ++it)
        std::cout << it->first << ",";
    std::cout << "\ncount:\t";
    for (auto it = freq_ctrs.begin(); it != freq_ctrs.end(); ++it)
        std::cout << it->second << ",";
    std::cout << std::endl;
}

}  // namespace priset
