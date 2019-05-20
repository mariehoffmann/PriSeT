#pragma once

#include <algorithm>
#include <fstream>
#include <string>

#include "../submodules/genmap/src/common.hpp"

namespace priset
{
/*
    using TLocations = std::map<seqan::Pair<TSeqNo, TSeqPos>,
             std::pair<std::vector<seqan::Pair<TSeqNo, TSeqPos> >,
                       std::vector<seqan::Pair<TSeqNo, TSeqPos> > > >;
                       */

// TODO: globally or hierarchical?
// pre-filter and sequence fetch
// 1. filter candidates by number of occurences only independent of their chemical suitability
// 2. fetch sequence and check chemical constraints that need to hold for a single primer
template<typename io_cfg_type, typename primer_cfg_type, typename TLocations, typename TKmerLocations, typename TSeqSize, typename TDirectoryInformation>
void pre_frequency_filter(io_cfg_type & io_cfg, primer_cfg_type & primer_cfg, TLocations & locations,
    TKmerLocations & kmer_locations, TSeqSize min_occ, TDirectoryInformation & directoryInformation)
{
    // TODO: verify if locations are sorted by SeqID
    using key_type = typename TLocations::key_type;
    //using TSeqNo = typename seqan::Value<key_type, 1>::Type;
    //using TSeqPos = typename seqan::Value<key_type, 2>::Type;
    std::set<key_type> seen;
    key_type current;
    // kmer locations: std::vector<std::pair<TSeq, std::vector<seqan::Pair<TSeqNo, TSeqPos> > > >
    // uniqueness indirectly preserved by (SeqNo, SeqPos) if list sorted lexicographically

    //TSeqSize chromosomeCount = 0;
    std::vector<std::pair<std::string, TSeqSize> > fastaFiles; // fasta file, cumulative nbr. of chromosomes

    std::string lastFastaFile = std::get<0>(retrieveDirectoryInformationLine(directoryInformation[0]));

    std::cout << "print dirInfo: " << std::endl;
    for (auto const & row : directoryInformation)
    { // TODO: continue here
        auto const line = retrieveDirectoryInformationLine(row);
        std::cout << std::get<0>(line) << ", " << std::get<1>(line) << std::endl;
        /*
        if (lastFastaFile != std::get<0>(line))
        {
            fastaFiles.push_back({lastFastaFile, chromosomeCount - 1});
            lastFastaFile = std::get<0>(line);
        }
        ++chromosomeCount;
        */
    }
    /*
    for (auto const & location : locations)
    {
        // TODO: check vector lengths
        auto const & kmerPos = location.first;
        auto const & plusStrandLoc = location.second.first;
        auto const & minusStrandLoc = location.second.second;

        std::cout << kmerPos.i1 << ',' << kmerPos.i2 << std::endl;
        TSeqSize i = 0;
        TSeqSize nbrChromosomesInPreviousFastas = 0;
        for (auto const & fastaFile : fastaFiles)
        {
            csvFile << ';';
            bool subsequentIterations = false;
            TSeqSize num_occ = 0;
            while (i < plusStrandLoc.size() && plusStrandLoc[i].i1 <= fastaFile.second)
            {
                if (subsequentIterations)
                    csvFile << '|'; // separator for multiple locations in one column
                seqan::Pair<TSeqSize, TSeqSize> loc{plusStrandLoc[i].i1 - nbrChromosomesInPreviousFastas, plusStrandLoc[i].i2};
                //csvFile << (plusStrandLoc[i].i1 - nbrChromosomesInPreviousFastas) << ',' << plusStrandLoc[i].i2;
                subsequentIterations = true;
                ++num_occ; // TODO: checkable without loop?
                ++i;
            }
            nbrChromosomesInPreviousFastas = fastaFile.second + 1;
        }

        if ((it->second).first.size() < min_occ)
            continue;
        current = it->first;
        if (seen.find(current) != seen.end())
            std::cout << "Location: (" << seqan::getValueI1(current) << ", " << seqan::getValueI2(current) << ") already seen.\n";
        seen.insert(current);
    }
    */
}

// post-filter candidates fulfilling chemical constraints by their relative frequency
template<typename TKmerLocations, typename TSeqNo>
void post_frequency_filter(TKmerLocations kmer_locations, TSeqNo occurrence_freq)
{

}

// filter k-mers by frequency and chemical properties
template<typename io_cfg_type, typename primer_cfg_type, typename TLocations, typename TKmerLocations, typename TDirectoryInformation>
void filter(io_cfg_type & io_cfg, primer_cfg_type & primer_cfg, TLocations & locations, TKmerLocations & kmer_locations, TDirectoryInformation & directoryInformation)
{
    using TSeqNo = typename seqan::Value<typename TLocations::key_type, 1>::Type;
    // get number of taxids with at least one accession assigned to it
    std::string acc_file_name = io_cfg.get_acc_file();
    std::string line;
    TSeqNo min_occ = 0;
    std::ifstream in;
    in.open(acc_file_name);
    while(!in.eof()) {
	       getline(in, line);
	       ++min_occ;
    }
    min_occ = (min_occ) ? min_occ - 1 : 0;
    in.close();
    std::cout << "number of lines " << min_occ;
    // scale to be lower frequency bound for filters
    min_occ = TSeqNo(float(min_occ) * primer_cfg.get_occurence_freq());
    std::cout << "lower freq bound = " << min_occ << std::endl;
    pre_frequency_filter<io_cfg_type, primer_cfg_type, TLocations, TKmerLocations, TSeqNo>(io_cfg, primer_cfg, locations, kmer_locations, min_occ, directoryInformation);

    post_frequency_filter(kmer_locations, primer_cfg.get_occurence_freq());

}

}  // namespace priset
