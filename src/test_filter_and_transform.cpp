#include <cassert>
#include <chrono>
#include <iostream>
#include <experimental/filesystem>
#include <fstream>
#include <regex>
#include <sys/wait.h>
#include <unistd.h>
#include <unordered_map>
#include <vector>

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

#include "argument_parser.hpp"
#include "filter.hpp"
#include "gui.hpp"
#include "primer_cfg_type.hpp"
#include "types.hpp"
#include "utilities.hpp"

namespace fs = std::experimental::filesystem;
using namespace priset;

// g++ ../PriSeT/src/test_filter_and_transform.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -lstdc++fs -DNDEBUG -O3 -I/Users/troja/include -L/Users/troja/lib -lsdsl -ldivsufsort -o test_filter_and_transform

struct setup
{
    // TODO: make this runnable with arbitrarily located build folders
    std::string lib_dir = (fs::canonical("../PriSeT/src/tests/library/one_seq")).string();
    std::string work_dir = (fs::canonical("../PriSeT/src/tests/work/one_seq")).string();
    priset::io_cfg_type io_cfg{};
    priset::primer_cfg_type primer_cfg{};
    priset::TKLocations kLocations{};

    setup()
    {
        // basic init
        int argc = 5;
        char * argv[5] = {"priset", "-l", &lib_dir[0], "-w", &work_dir[0]};
        for (auto i = 0; i < argc; ++i)
            std::cout << "arg[" << i << "] = " << argv[i] << std::endl;
        priset::options opt(argc, argv, primer_cfg, io_cfg);
    }
};

/*
 * Passes chemical filter for Tm, CG content and runs for length 23, 21, 19
 * s1 = "AACGTAACGTAACGTACGTACGT" at pos 10
 *       01234567890123456789012
 * k1_pattern_2 = 000101010000|encode(s1)
 *
 */
void test_one_seq()
{
    setup su{};
    priset::TKLocations locations;
    std::vector<TLocation> loc{TLocation{1, 10}}, vd{}; // seqan::Pair<priset::TSeqNo, priset::TSeqPos>
    using TLocationPair = std::pair<std::vector<TLocation>, std::vector<TLocation>>;
    locations.insert({TKLocation{1, 10, 19}, TLocationPair{loc, vd}});
    locations.insert({TKLocation{1, 10, 21}, TLocationPair{loc, vd}});
    locations.insert({TKLocation{1, 10, 23}, TLocationPair{loc, vd}});

    priset::TReferences references;
    priset::TKmerIDs kmerIDs;
    priset::TSeqNoMap seqNoMap;
    priset::TSeqNo cutoff = 1;
    TKmerCountStats kmer_count_stats;
    filter_and_transform(su.io_cfg, su.primer_cfg, locations, references, kmerIDs, seqNoMap, cutoff, kmer_count_stats);
    std::cout << "References:\n";
    for (TReference reference : references)
    {
        std::cout << "\t";
        for (auto bit : reference)
            std::cout << bit;
        std::cout << std::endl;
    }
    std::cout << "KmerIDs:\n";
    for (auto kmerID_list : kmerIDs)
    {
        for (TKmerID kmerID : kmerID_list)
            std::cout << kmerID << " | ";
        std::endl;
    }
}

int main()
{

    test_one_seq();
}
