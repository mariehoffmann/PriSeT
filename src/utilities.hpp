#pragma once

#include <iostream>

#include <seqan/basic.h>

#include "types.hpp"

/*
 * TLocations = std::map<seqan::Pair<TSeqNo, TSeqPos>,
 *         std::pair<   std::vector<seqan::Pair<TSeqNo, TSeqPos> >,
 *                      std::vector<seqan::Pair<TSeqNo, TSeqPos> > > >
 */
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

// Retrieve DNA sequence from txt.concat given a set of locations
template<typename primer_cfg_type>
void lookup_sequences(priset::TKmerLocations & kmer_locations, priset::io_cfg_type & io_cfg, primer_cfg_type & primer_cfg_type)
{
    std::string fastaFile = std::get<0>(retrieveDirectoryInformationLine(directoryInformation[0]));
    TSeqNo startPos = 0;
    TSeqNo seqLen = 0;
    TSeqPos offset;
    TKmerLocations::const_iterator it = kmer_locations.begin();
    TSeqNo next_kmer_id = seqan::getValueI1<TSeqNo, TSeqPos>(it->second[0]);

    TSeqPos kmer_length = primer_cfg.get_kmer_length();
    for (uint64_t id = 0; id < length(directoryInformation); ++id)
    {
        auto const row = retrieveDirectoryInformationLine(directoryInformation[id]);
        seqLen = std::get<1>(row);
        auto chromosomeNames = std::get<2>(row);

        // std::vector<std::pair<TSeq, std::vector<seqan::Pair<priset::TSeqNo, priset::TSeqPos> > > >
        // if sequence id corresponds to first sequence number, we grep the k-mer sequence
        if (id == next_kmer_id)
        {
            offset = seqan::getValueI2<TSeqNo, TSeqPos>(it->second[0]);
            auto const & kmer = infixWithLength(text.concat, startPos + offset, startPos + offset + seqLen);
            // copy kmer into first position of locations vector
            (*it).first = kmer;
            // forward next kmer iterator
            ++it;
            if (it == kmer_locations.end()) // no more kmers to resolve
                break
            next_kmer_id = seqan::getValueI1<TSeqNo, TSeqPos>(it->second[0]);
        }
        if (std::get<0>(row) != fastaFile)
            break;
        startPos += seqLen;
    }
}

// set directory information as needed by genmap's fasta file parsing
template<typename TDirectoryInformation>
void set_directoryInformation(std::string & index_path_base_ids, TDirectoryInformation & directoryInformation)
{
    seqan::open(directoryInformation, seqan::toCString(index_path_base_ids), seqan::OPEN_RDONLY);
    // dummy entry enforces that the mappability is computed for the last file in the while loop.
    seqan::appendValue(directoryInformation, "dummy.entry;0;chromosomename");
}
