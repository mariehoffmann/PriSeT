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

#include "errors.hpp"
#include "filter.hpp"
#include "fm.hpp"
#include "gui.hpp"
#include "io_cfg_type.hpp"
#include "primer_cfg_type.hpp"
#include "taxonomy.hpp"
#include "types.hpp"
#include "utilities.hpp"

/*
 * usage        g++ ../PriSeT/src/priset.cpp -Wno-write-strings -std=c++17 -lstdc++fs -Wall -Wextra -o priset
 *              ./priset <lib_dir> <work_dir>
 * e.g.         ./priset ~/priset/library ~/priset/work
 *
 * src_dir      path to folder containing fasta (*.fa) and taxonomy file (*.tax)
 * work_dir     path to store indices, mappings, annotations, and other results
 */
int main(int argc, char** argv)
{
    if (argc < 3)
        std::cout << "ERROR: " << ARG_ERROR << std::endl, exit(0);

    // set path prefixes for library files
    priset::io_cfg_type io_cfg{argv[1], argv[2]};

    // get instance to primer sequence settings
    priset::primer_cfg_type primer_cfg{};

    // build taxonomy in RAM
    priset::taxonomy tax{io_cfg.get_tax_file()};
    //tax.print_taxonomy();

    // create FM index
    priset::fm_index<priset::io_cfg_type>(io_cfg);

    // dictionary for storing FM mapping results
    priset::TLocations locations;

    // directory info needed for genmap's fasta file parser
    priset::TDirectoryInformation directoryInformation;

    // container for fasta header lines
    priset::TSequenceNames sequenceNames;

    // container for fasta sequence lengths
    priset::TSequenceLengths sequenceLengths;

    // compute k-mer mappings
    priset::fm_map2(io_cfg, primer_cfg, locations, directoryInformation); //, sequenceNames, sequenceLengths);
    priset::print_locations(locations);

    // filter k-mers by frequency and chemical properties
    // TODO: result structure for references and k-mer pairs: candidates/matches
    // vector storing k-mer IDs and their locations, i.e. {TSeq: [(TSeqAccession, TSeqPos)]}
    priset::TKmerLocations kmer_locations;
    // dictionary to resolve kmer IDs and their sequences
    priset::TKmerMap kmer_map;
    priset::pre_filter_main<priset::TSequenceNames, priset::TSequenceLengths>(io_cfg, primer_cfg, locations, kmer_locations, kmer_map, directoryInformation); //, sequenceNames, sequenceLengths);
    // TODO: delete locations
    priset::TKmerPairs pairs;
    priset::combine(primer_cfg, kmer_locations, kmer_map, pairs);
    // test chemical constraints of pairs and filter
    //priset::post_filter_main(primer_cfg, kmer_locations, pairs);
    priset::create_table(io_cfg, kmer_locations, pairs);
    // create app script
    if (! priset::gui::generate_app(io_cfg) && priset::gui::compile_app(io_cfg))
        std::cout << "ERROR: gui::generate_app or gui::compile_app returned false\n";

    return 0;
}
