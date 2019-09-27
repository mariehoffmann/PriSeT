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

#include "../argument_parser.hpp"
#include "../combine_types.hpp"
#include "../filter.hpp"
#include "../gui.hpp"
#include "../primer_cfg_type.hpp"
#include "../types.hpp"
#include "../utilities.hpp"

namespace fs = std::experimental::filesystem;
using namespace priset;

// g++ ../PriSeT/src/tests/filter_and_transform_test.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -lstdc++fs -DNDEBUG -O3 -I/Users/troja/include -L/Users/troja/lib -lsdsl -ldivsufsort -o filter_and_transform_test

struct setup
{
    // TODO: make this runnable with arbitrarily located build folders
    std::string lib_dir = (fs::canonical("../PriSeT/src/tests/library/one_seq")).string();
    std::string work_dir = (fs::canonical("../PriSeT/src/tests/work/one_seq")).string();
    io_cfg_type io_cfg{};
    primer_cfg_type primer_cfg{};
    TKLocations locations{};
    using TLocationPair = std::pair<std::vector<TLocation>, std::vector<TLocation>>;

    TReferences references;
    TKmerIDs kmerIDs;
    TSeqNoMap seqNoMap;
    TSeqNo cutoff = 1;
    TKmerCounts kmerCounts;

    fs::path idx_dir = work_dir + "/index";
    fs::path idx_zip = work_dir + "/index.zip";
    fs::path tmp_dir = work_dir + "/tmp";

    setup()
    {
        // delete previous index dir if existent
        fs::remove_all(idx_dir);
        // unzip index.zip into same named directory
        std::system(("unzip -n -d " + work_dir + " " + idx_zip.string()).c_str());

        // create tmp dir
        if (fs::exists(tmp_dir))
        {
            std::cout << "tmp_dir exists, delete it\n";
            fs::remove_all(tmp_dir);
        }
        if (!fs::create_directory(tmp_dir))
            std::cout << "ERROR: could not create tmp_dir = " << tmp_dir << std::endl;
        std::cout << "tmp_dir in setup = " << tmp_dir << std::endl;

        std::cout << "lib_dir in setup = " << lib_dir << std::endl;
        std::cout << "work_dir in setup = " << work_dir << std::endl;

        int argc = 6;
        char * argv[6] = {"priset", "-l", &lib_dir[0], "-w", &work_dir[0], "-s"};
        for (auto i = 0; i < argc; ++i)
            std::cout << "arg[" << i << "] = " << argv[i] << std::endl;
        options opt(argc, argv, primer_cfg, io_cfg);

        std::vector<TLocation> loc1{TLocation{0, 10}}, vd{}; // seqan::Pair<TSeqNo, TSeqPos>
        locations.insert({TKLocation{0, 10, 19}, TLocationPair{loc1, vd}});
        locations.insert({TKLocation{0, 10, 21}, TLocationPair{loc1, vd}});
        locations.insert({TKLocation{0, 10, 23}, TLocationPair{loc1, vd}});

        std::vector<TLocation> loc2{TLocation{0, 90}}; // seqan::Pair<TSeqNo, TSeqPos>
        locations.insert({TKLocation{0, 90, 21}, TLocationPair{loc2, vd}});

    }
};

/*
 * tests interfaces between fm_map and filter_and_transform
 * Passes chemical filter for Tm, CG content and runs for length 23, 21, 19
 *
 * kmer 1-3 at position 10: "(C)AACGTAACGTAACGTACGTACGT"
 * head = (1 << (63-3)) + (1 << (63-5)) + (1 << (63-7)) = 288238240355526424
 *
 * kmer 2 at position 90: "TAGCTAACTACATAGCTACGA"
 * head = (1 << (63 - (21-16)))                          = 288238240355526424
 *
    dTm(1513281700780251931, 19, 288238240355526424, 21) = 12

    AACGTAACGTAACGTACGT         AT_cnt = 11 CG_cnt = 8 => 11*2 + 8*4 = 54 degrees
    TAGCTAACTACATAGCTACGA       AT_cnt = 13 CG_cnt = 8 => 26 + 32 = 58 degrees
 */
void test_filter_and_transform()
{
    setup su{};

    filter_and_transform(su.io_cfg, su.primer_cfg, su.locations, su.references, su.kmerIDs, su.seqNoMap, su.cutoff, su.kmerCounts);
    std::cout << "References:\n";
    for (TReference reference : su.references)
    {
        std::cout << "\t";
        for (auto bit : reference)
            std::cout << bit;
        std::cout << std::endl;
    }
    std::cout << "KmerIDs:\n";
    for (auto kmerID_list : su.kmerIDs)
    {
        for (TKmerID kmerID : kmerID_list)
            std::cout << kmerID << " | ";
        std::cout << std::endl;
    }

    TKmerID expect1 = dna_encoder("AACGTAACGTAACGTACGTACGT") + (ONE_LSHIFT_63 >> 3) + (ONE_LSHIFT_63 >> 5) + (ONE_LSHIFT_63 >> 7);
    TKmerID expect2 = dna_encoder("TAGCTAACTACATAGCTACGA") + (ONE_LSHIFT_63 >> 5);
    if (!su.kmerIDs.size() || su.kmerIDs[0].size() != 2 || su.kmerIDs[0][0] != expect1 || su.kmerIDs[0][1] != expect2)
        std::cout << "ERROR: expect kmerID1 = " << expect1 << " and kmerID2 = " << expect2 << ", but got nothing or a wrong kmerID\n";
    else
        std::cout << "SUCCESS: Result as expected!\n";
}

void test_combine()
{
    setup su{};
    filter_and_transform(su.io_cfg, su.primer_cfg, su.locations, su.references, su.kmerIDs, su.seqNoMap, su.cutoff, su.kmerCounts);
    //using TPair = TPair;
    using TPairList = TPairList<TPair<TCombinePattern<TKmerID, TKmerLength>>>;
    TPairList pairs;
    //std::vector<_Ch_type, std::allocator<_CharT> > >(priset::primer_cfg_type&, priset::TKmerIDs&, priset::TPairList<priset::TPair<priset::TCombinePattern<long long unsigned int, long long int> > >&)'
//     print_combinations<>(su.primer_cfg, su.kmerIDs, pairs);
    combine(su.primer_cfg, su.references, su.kmerIDs, pairs, su.kmerCounts);
    print_combinations<TPairList>(su.kmerIDs, pairs);
}

int main()
{
    test_filter_and_transform();
    //test_combine();
    return 0;
}
