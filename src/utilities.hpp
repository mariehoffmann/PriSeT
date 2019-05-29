#pragma once

#include <iostream>

#include <seqan/basic.h>

#include "chemistry.hpp"
#include "types.hpp"

namespace priset
{
template<typename TLocations>
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

// forward declaration
struct primer_cfg_type;

// Retrieve DNA sequence from txt.concat given a set of locations
// lookup_sequences<primer_cfg_type>(kmer_locations, io_cfg, primer_cfg, directoryInformation);
template<typename primer_cfg_type>
void lookup_sequences(TKmerLocations & kmer_locations, io_cfg_type const & io_cfg, primer_cfg_type const & primer_cfg, TDirectoryInformation const & directoryInformation)
{
    std::string fastaFile = std::get<0>(retrieveDirectoryInformationLine(directoryInformation[0]));
    TSeqNo startPos = 0;    // accession separation offset
    TSeqPos offset;         // kmer offset within accession
    TSeqPos kmer_length = primer_cfg.get_kmer_length();
    TKmerLocations::iterator kmer_it = kmer_locations.begin();
    // unique identifier to be assigned to TKmer struct
    TKmerID ID = 0;
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
        // kmer is struct
        (*kmer_it).first = TKmer{str, ID++, chemistry::get_Tm(primer_cfg, str)}; // direct assignment of kmer possible?
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

}  // namespace priset
