#include <array>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <experimental/filesystem>
#include <fstream>
#include <numeric>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

#include "../src/argument_parser.hpp"
#include "../src/combine_types.hpp"
#include "../src/filter.hpp"
#include "../src/fm.hpp"
#include "../src/io_cfg_type.hpp"
#include "../src/primer_cfg_type.hpp"
#include "../src/priset.hpp"
#include "../src/types.hpp"
#include "../src/utilities.hpp"

namespace fs = std::experimental::filesystem;

// g++ ../PriSeT/tests/clade_3041.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -lstdc++fs -Wno-unknown-pragmas -lstdc++fs -DNDEBUG -O3 -I/Users/troja/include -L/Users/troja/lib -lsdsl -ldivsufsort -o clade_3041

struct setup
{
    std::string lib_dir; // = "tactac/subset/3041").string();
    //std::cout << "lib_dir = " << lib_dir << std::endl;
    std::string work_dir; // = fs::canonical("/Volumes/plastic_data/priset//work/3041").string();

    fs::path idx_dir;
    fs::path idx_zip;
    fs::path tmp_dir;

    setup(std::string lib_dir, std::string work_dir)
    {
        lib_dir = fs::canonical(lib_dir).string();
        work_dir = fs::canonical(work_dir).string();
        idx_dir = work_dir + "/index";
        tmp_dir = work_dir + "/tmp";

        std::cout << "lib_dir in setup = " << lib_dir << std::endl;
        std::cout << "work_dir in setup = " << work_dir << std::endl;
    }
};

int main(int argc, char ** argv)
{
    if (argc != 3)
    {
        std::cout << "Give paths to lib and work dirs.\n";
        exit(-1);
    }
    setup su{argv[1], argv[2]};

    unsigned const priset_argc = 6;
    char * const priset_argv[priset_argc] = {"priset", "-l", argv[1], "-w", argv[2], "-s"};
    for (unsigned i = 0; i < priset_argc; ++i) std::cout << priset_argv[i] << " ";
    std::cout << std::endl;

    // collect number of kmers or kmer pairs left after relevant processing steps
    TKmerCounts kmerCounts{0, 0, 0, 0};

    // set path prefixes for library files
    io_cfg_type io_cfg{};

    // get instance to primer sequence settings
    primer_cfg_type primer_cfg{};

    // parse options and init io and primer configurators
    options opt(priset_argc, priset_argv, primer_cfg, io_cfg);

    std::chrono::time_point<std::chrono::system_clock> start, finish;

    TKLocations locations;
    TDirectoryInformation directoryInformation;
    TSequenceNames sequenceNames;
    TSequenceLengths sequenceLengths;

    // compute k-mer mappings
    start = std::chrono::high_resolution_clock::now();
    fm_map(io_cfg, primer_cfg, locations);
    finish = std::chrono::high_resolution_clock::now();
    std::cout << "Runtime fm_map: " << std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();

    std::cout << "INFO: kmers init = " << locations.size() << std::endl;

    /*for (auto it = locations.begin(); it != locations.end(); ++it)
    {
        std::cout << "(" << std::get<0>(it->first) << ", " << std::get<1>(it->first) << ", K = " << std::get<2>(it->first) << "): [";
        if (it->second.first.size() > 1)
            std::cout << "(" << it->second.first[0].i1 << ", " << it->second.first[0].i2 << "),(" << it->second.first[1].i1 << ", " << it->second.first[0].i2 << ")\n";
    }*/
    TReferences references;
    TKmerIDs kmerIDs;
    TSeqNoMap seqNoMap;
    start = std::chrono::high_resolution_clock::now();
    filter_and_transform(io_cfg, locations, references, kmerIDs, seqNoMap, kmerCounts);
    finish = std::chrono::high_resolution_clock::now();
    std::cout << "Runtime filter_and_transform: " << std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count()  << std::endl;

    std::cout << "INFO: kmers after filter1 & transform = " << get_num_kmers(kmerIDs) << std::endl;

    // TODO: delete locations
    using TPairList = TPairList<TPair<TCombinePattern<TKmerID, TKmerLength>>>;
    TPairList pairs;

    start = std::chrono::high_resolution_clock::now();
    combine<TPairList>(references, kmerIDs, pairs, kmerCounts);
    finish = std::chrono::high_resolution_clock::now();
    std::cout << "Runtime combine: " << std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count()  << std::endl;

    std::cout << "INFO: pairs after combiner = " << get_num_pairs<TPairList>(pairs) << std::endl;

    start = std::chrono::high_resolution_clock::now();
    filter_pairs(references, kmerIDs, pairs);
    finish = std::chrono::high_resolution_clock::now();
    std::cout << "Runtime filter_pairs: " << std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count() << std::endl;

    std::cout << "INFO: pairs after frequency cutoff = " << get_num_pairs<TPairList>(pairs) << std::endl;

    // collect unique primer sequences
    std::unordered_set<uint64_t> kmers_unique;

    for (TPair<TCombinePattern<TKmerID, TKmerLength>> pair : pairs)
    {
        for (uint8_t i = 0; i < 100; ++i)
        {
            TKmerID kmer_fwd = kmerIDs[pair.reference][pair.r_fwd];
            TKmerID kmer_rev = kmerIDs[pair.reference][pair.r_rev];

            if (pair.cp[i])
            {
                kmers_unique.insert(get_code(kmer_fwd, ONE_LSHIFT_63 >> (i/10)));
                kmers_unique.insert(get_code(kmer_rev, ONE_LSHIFT_63 >> (i % 10)));
            }
        }
    }
    // write primers into file
    fs::path primer_file = su.tmp_dir / "primers_3041.csv";

    std::ofstream ofs;
    ofs.open(primer_file);
    ofs << "primer\n";
    for (auto primer : kmers_unique)
        ofs << dna_decoder(primer) << "\n";
    ofs.close();
    return 0;
}
