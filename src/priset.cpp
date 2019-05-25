// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================
//          Author: Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
//          Manual: https://github.com/mariehoffmann/PriSeT

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <iostream>
#include <unistd.h>
#include <vector>

#include <seqan/basic.h>

#include "display.hpp"
#include "errors.hpp"
#include "filter.hpp"
#include "fm.hpp"
#include "io_config.hpp"
#include "primer_config.hpp"
#include "taxonomy.hpp"
#include "types.hpp"
#include "utilities.hpp"

/*
 * usage        g++ priset.cpp -Wno-write-strings -std=c++17 -lstdc++fs -Wall -Wextra -o priset
 *              ./priset <lib_dir> <work_dir>
 * e.g.         ./src/priset ~/priset/library ~/priset/work
 *
 * src_dir      path to folder containing fasta (*.fa) and taxonomy file (*.tax)
 * work_dir     path to store indices, mappings, annotations, and other results
 */
int main(int argc, char** argv)
{
    if (argc < 3)
        std::cout << "ERROR: " << ARG_ERROR << std::endl, exit(0);

    // set path prefixes for library files
    priset::io_config io_cfg{argv[1], argv[2]};

    // get instance to primer sequence settings
    priset::primer_config<priset::TSeq> primer_cfg{};

    // build taxonomy in RAM
    priset::taxonomy tax{io_cfg.get_tax_file()};
    //tax.print_taxonomy();

    // create FM index
    priset::fm_index<priset::io_config>(io_cfg);

    // dictionary for storing FM mapping results
    priset::TLocations locations;

    // directory info needed for genmap's fasta file parser
    priset::TDirectoryInformation directoryInformation;

    // container for fasta header lines
    priset::TSequenceNames sequenceNames;

    // container for fasta sequence lengths
    priset::TSequenceLengths sequenceLengths;

    // compute k-mer mappings
    priset::fm_map<priset::io_config, priset::primer_config<priset::TSeq>, priset::TLocations, priset::TDirectoryInformation, priset::TSequenceNames, priset::TSequenceLengths>(io_cfg, primer_cfg, locations, directoryInformation, sequenceNames, sequenceLengths);
    print_locations(locations);

    // filter k-mers by frequency and chemical properties
    // TODO: result structure for references and k-mer pairs: candidates/matches
    // dictionary for storing k-mers and their locations, i.e. {TSeq: [(TSeqAccession, TSeqPos)]}
    priset::TKmerLocations kmer_locations;
    priset::filter<priset::primer_config<priset::TSeq>, priset::TLocations, priset::TKmerLocations, priset::TDirectoryInformation, priset::TSequenceNames, priset::TSequenceLengths>(io_cfg, primer_cfg, locations, kmer_locations, directoryInformation, sequenceNames, sequenceLengths);
    // TODO: delete locations

    // display
    priset::display(io_cfg/*, candidates*/);

    return 0;
}
