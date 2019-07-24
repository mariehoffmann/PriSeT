#pragma once

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <iostream>
#include <unistd.h>
#include <vector>

#include <seqan/basic.h>

#include "argument_parser.hpp"
#include "errors.hpp"
#include "filter.hpp"
#include "fm.hpp"
#include "gui.hpp"
#include "io_cfg_type.hpp"
#include "primer_cfg_type.hpp"
#include "taxonomy.hpp"
#include "types.hpp"
#include "utilities.hpp"

#define DEBUG 1

/*
 * usage        g++ ../PriSeT/src/priset.cpp -Wno-write-strings -std=c++17 -lstdc++fs -Wall -Wextra -o priset
 *              ./priset <lib_dir> <work_dir>
 * e.g.         ./priset -l ~/priset/library -w ~/priset/work
 *
 * src_dir      path to folder containing fasta (*.fa) and taxonomy file (*.tax)
 * work_dir     path to store indices, mappings, annotations, and other results
 */
int priset_main(int argc, char * const * argv, std::array<size_t, priset::TIMEIT::SIZE> * runtimes = nullptr)
{
    bool timeit_flag(runtimes);
    // Store start and finish times for optional runtime measurements
    std::chrono::time_point<std::chrono::system_clock> start, finish;
    if (timeit_flag)  // init timeit array
        for (auto i = 0; i < priset::TIMEIT::SIZE; ++i)
            runtimes->at(i) = 0;

    // set path prefixes for library files
    priset::io_cfg_type io_cfg{};

    // get instance to primer sequence settings
    priset::primer_cfg_type primer_cfg{};

    // parse options and init io and primer configurators
    priset::options opt(argc, argv, primer_cfg, io_cfg);

    // build taxonomy in RAM
    //priset::taxonomy tax{io_cfg.get_tax_file()};
    //tax.print_taxonomy();

    // create FM index if SKIP_IDX not in argument list
    int ret_code;
    if (io_cfg.skip_idx())
    {
        std::cout << "MESSAGE: skip index recomputation" << std::endl;
    }
    else if ((ret_code = priset::fm_index(io_cfg)))
    {
        std::cout << "ERROR: " << ret_code << std::endl;
        exit(-1);
    }
    // quit here for index computation without subsequent mappability
    if (io_cfg.idx_only())
    {
        std::cout << "MESSAGE: index recomputation only" << std::endl;
        return 0;
    }

    // dictionary for storing FM mapping results
    priset::TLocations locations;

    // directory info needed for genmap's fasta file parser
    priset::TDirectoryInformation directoryInformation;

    // container for fasta header lines
    priset::TSequenceNames sequenceNames;

    // container for fasta sequence lengths
    priset::TSequenceLengths sequenceLengths;

    if (timeit_flag)
        start = std::chrono::high_resolution_clock::now();
    // compute k-mer mappings
    priset::fm_map(io_cfg, primer_cfg, locations);
    if (timeit_flag)
    {
        finish = std::chrono::high_resolution_clock::now();
        // duration obj
        //auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(finish - start);
        runtimes->at(priset::TIMEIT::MAP) = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
    }

    // Do not modify or delete STATS lines, since they are captured for statistical analysis
    std::cout << "INFO: kmers init = " << locations.size() << std::endl;

    // filter k-mers by frequency and chemical properties
    // TODO: result structure for references and k-mer pairs: candidates/matches
    // vector storing k-mer IDs and their locations, i.e. {TSeq: [(TSeqAccession, TSeqPos)]}
    priset::TKmerLocations kmer_locations;
    // dictionary to resolve kmer IDs and their sequences
    priset::TKmerMap kmer_map;
    if (timeit_flag)
        start = std::chrono::high_resolution_clock::now();
    priset::pre_filter_main(io_cfg, primer_cfg, locations, kmer_locations, kmer_map);
    if (timeit_flag)
    {
        finish = std::chrono::high_resolution_clock::now();
        runtimes->at(priset::TIMEIT::FILTER1) = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
    }

    // TODO: delete locations
    priset::TKmerPairs pairs;
    if (timeit_flag)
        start = std::chrono::high_resolution_clock::now();
    priset::combine(primer_cfg, kmer_locations, kmer_map, pairs);
    if (timeit_flag)
    {
        finish = std::chrono::high_resolution_clock::now();
        runtimes->at(priset::TIMEIT::COMBINER) = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
    }

    std::cout << "INFO: pairs combined = " << pairs.size() << std::endl;

    // test chemical constraints of pairs and filter
    if (timeit_flag)
        start = std::chrono::high_resolution_clock::now();

    //priset::post_filter_main(primer_cfg, kmer_locations, pairs);
    if (timeit_flag)
    {
        finish = std::chrono::high_resolution_clock::now();
        runtimes->at(priset::TIMEIT::FILTER2) = std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
    }

    std::cout << "INFO: pairs filtered = " << pairs.size() << std::endl;
    if (!timeit_flag)
    {
        priset::create_table(io_cfg, kmer_locations, pairs, kmer_map);
        // create app script
        if (! priset::gui::generate_app(io_cfg) && priset::gui::compile_app(io_cfg))
            std::cout << "ERROR: gui::generate_app or gui::compile_app returned false\n";
    }
    return 0;
}
