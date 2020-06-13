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

// macro used in order to ignore frequency cutoff filter for this test
#define PRISET_TEST

#include "../src/argument_parser.hpp"
#include "../src/filter.hpp"
#include "../src/gui.hpp"
#include "../src/PrimerConfig.hpp"
#include "../src/types.hpp"
#include "../src/utilities.hpp"

namespace fs = std::experimental::filesystem;
using namespace priset;


// g++ ../PriSeT/tests/transform_and_filter_test.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -lstdc++fs -DNDEBUG -O3 -I/Users/troja/include -L/Users/troja/lib -lsdsl -ldivsufsort -o transform_and_filter_test

struct setup
{
    // TODO: make this runnable with arbitrarily located build folders
    std::string lib_dir = (fs::canonical("../PriSeT/tests/library/one_seq")).string();
    std::string work_dir = (fs::canonical("../PriSeT/tests/work/one_seq")).string();
    IOConfig io_cfg{};
    PrimerConfig primer_cfg{};
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

void test_transform_and_filter()
{
    setup su{};
    transform_and_filter(su.io_cfg, su.primer_cfg, su.locations, su.references, su.kmerIDs, su.seqNoMap, su.kmerCounts);
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

    // k = 23: AT = 13, CG = 10, Tm = 26 + 40 = 66!, CG_content = .43, no CG clamp
    // k = 21: AT = 12, CG = 9, Tm = 24 + 36 = 60, CG_content = .42, no CG clamp
    // k = 19: AT = 11, CG = 8, Tm = 22 + 32 = 54, CG_content = .42, no CG clamp
    TKmerID expect1 = dna_encoder("AACGTAACGTAACGTACGTAC") + (ONE_LSHIFT_63 >> 3) + (ONE_LSHIFT_63 >> 5);
    // k = 21: AT = 12, CG = 9, Tm = 24 + 36 = 60, CG_content = .43, no CG clamp
    TKmerID expect2 = dna_encoder("TAGCTAACTACATAGCTACGC") + (ONE_LSHIFT_63 >> 5);

    if (!su.kmerIDs.size() || su.kmerIDs[0].size() != 2 || su.kmerIDs[0][0] != expect1 || su.kmerIDs[0][1] != expect2)
        std::cout << "ERROR: expect kmerID1 = " << expect1 << " and kmerID2 = " << expect2 << ", but got nothing or a wrong kmerID\n";
    else
        std::cout << "SUCCESS: Result as expected!\n";
}

void test_combine()
{
    setup su{};
    transform_and_filter(su.io_cfg, su.primer_cfg, su.locations, su.references, su.kmerIDs, su.seqNoMap, su.kmerCounts);
    //using TPair = TPair;
    using TPairList = TPairList<TPair<CombinePattern>>;
    TPairList pairs;
    //std::vector<_Ch_type, std::allocator<_CharT> > >(priset::PrimerConfig&, priset::TKmerIDs&, priset::TPairList<priset::TPair<priset::CombinePattern>> &)'
//     print_combinations<>(su.primer_cfg, su.kmerIDs, pairs);
    combine(su.references, su.kmerIDs, pairs, su.kmerCounts);
    print_combinations<TPairList>(su.kmerIDs, pairs);
}

int main()
{
    test_transform_and_filter();
    //test_combine();
    return 0;
}
