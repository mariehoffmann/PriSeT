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
    io_cfg_type io_cfg{};
    primer_cfg_type primer_cfg{};
    TKLocations kLocations{};

    fs::path idx_dir = work_dir + "/index";
    fs::path idx_zip = work_dir + "/index.zip";
    fs::path tmp_dir = work_dir + "/tmp";

    setup()
    {
        // unzip index.zip into same named directory
        std::system(("unzip -n -d " + work_dir + " " + idx_zip.string()).c_str());
        // create tmp dir
        if (!fs::create_directory(tmp_dir))
            std::cout << "ERROR: could not create tmp_dir = " << tmp_dir << std::endl;
        std::cout << "lib_dir in setup = " << lib_dir << std::endl;
        std::cout << "work_dir in setup = " << work_dir << std::endl;

        int argc = 6;
        char * argv[6] = {"priset", "-l", &lib_dir[0], "-w", &work_dir[0], "-s"};
        for (auto i = 0; i < argc; ++i)
            std::cout << "arg[" << i << "] = " << argv[i] << std::endl;
        options opt(argc, argv, primer_cfg, io_cfg);

    }
};

/*
 * tests interfaces between fm_map and filter_and_transform
 * Passes chemical filter for Tm, CG content and runs for length 23, 21, 19
 * s1 = "AACGTAACGTAACGTACGTACGT" at pos 10
 *       01234567890123456789012
 * k1_pattern_2 = 000101010000_{2}|133286677332880_{10} = 47421082764723088
 *
 */
void test_one_seq()
{
    setup su{};
    TKLocations locations; // zero-based
    std::vector<TLocation> loc{TLocation{0, 10}}, vd{}; // seqan::Pair<TSeqNo, TSeqPos>
    using TLocationPair = std::pair<std::vector<TLocation>, std::vector<TLocation>>;
    // encode("AACGTAACGTAACGTACGT") = 520648278928
    locations.insert({TKLocation{0, 10, 19}, TLocationPair{loc, vd}});
    // encode("AACGTAACGTAACGTACGTAC") = 5743328510864
    locations.insert({TKLocation{0, 10, 21}, TLocationPair{loc, vd}});
    // encode("AACGTAACGTAACGTACGTACGT") = 133286677332880
    locations.insert({TKLocation{0, 10, 23}, TLocationPair{loc, vd}});

    TReferences references;
    TKmerIDs kmerIDs;
    TSeqNoMap seqNoMap;
    TSeqNo cutoff = 1;
    TKmerCounts kmerCounts;
    filter_and_transform(su.io_cfg, su.primer_cfg, locations, references, kmerIDs, seqNoMap, cutoff, kmerCounts);
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
        std::cout << std::endl;
    }
}

int main()
{

    test_one_seq();
    return 0;
}
