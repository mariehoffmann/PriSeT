#include <array>
#include <chrono>
// #include <cstdlib>
#include <iostream>
#include <experimental/filesystem>
// #include <fstream>
// #include <numeric>
#include <sys/wait.h>
#include <unistd.h>
#include <unordered_map>
#include <vector>

#include <seqan/basic.h>

#include "../src/argument_parser.hpp"
#include "../src/errors.hpp"
#include "../src/filter.hpp"
#include "../src/fm.hpp"
#include "../src/gui.hpp"
#include "../src/io_cfg_type.hpp"
#include "../src/output.hpp"
#include "../src/primer_cfg_type.hpp"
#include "../src/taxonomy.hpp"
#include "../src/types.hpp"
#include "../src/utilities.hpp"

namespace fs = std::experimental::filesystem;
using namespace priset;

// static constexpr bool outputProgress = false;

// -mpopcnt gives speedup of 5, otherwise call __popcountdi2 called for __builtin_popcountll
// g++ ../PriSeT/tests/performance_test.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -lstdc++fs -Wno-unknown-pragmas -lstdc++fs -DNDEBUG -O3 -mpopcnt -I/Users/troja/include -L/Users/troja/lib -lsdsl -ldivsufsort -o performance_test

struct setup
{
    std::string lib_dir = fs::canonical("../PriSeT/tests/library/3041").string();
    //std::cout << "lib_dir = " << lib_dir << std::endl;
    std::string work_dir = fs::canonical("../PriSeT/tests/work/3041").string();

    fs::path idx_dir = work_dir + "/index";
    fs::path idx_zip = work_dir + "/index.zip";
    fs::path tmp_dir = work_dir + "/tmp";

    setup()
    {
        // unzip index.zip into same named directory
        std::system(("unzip -n -d " + work_dir + " " + idx_zip.string()).c_str());
        // clear and create tmp dir
        std::cout << "tmpdir exists: " << fs::exists(tmp_dir) << std::endl;
        if (fs::exists(tmp_dir))
        {
            std::cout << "tmp_dir exists, delete it\n";
            fs::remove_all(tmp_dir);
        }
        if (!fs::create_directory(tmp_dir))
            std::cout << "ERROR: could not create tmp_dir = " << tmp_dir << std::endl;
        std::cout << "lib_dir in setup = " << lib_dir << std::endl;
        std::cout << "work_dir in setup = " << work_dir << std::endl;

    }

    void down()
    {
        // delete index dir
        if (fs::remove_all(idx_dir))
            std::cout << "ERROR: could not remove idx_dir = " << idx_dir << std::endl;
        // delete tmp dir
        if (fs::remove_all(tmp_dir))
            std::cout << "ERROR: could not remove tmp_dir = " << tmp_dir << std::endl;
    }
};

/* Measure runtime for PriSeT components */
int main(/*int argc, char ** argv*/)
{
    setup su{};
    std::array<size_t, priset::TIMEIT::SIZE> runtimes;

    unsigned const argc = 6;
    char * const argv[argc] = {"priset", "-l", &su.lib_dir[0], "-w", &su.work_dir[0], "-s"};
    for (unsigned i = 0; i < argc; ++i) std::cout << argv[i] << " ";
    std::cout << std::endl;

    ///////////////////////////////////////////
    // Store start and finish times for optional runtime measurements
    std::chrono::time_point<std::chrono::system_clock> start, finish;

    // collect number of kmers or kmer pairs left after relevant processing steps
    TKmerCounts kmerCounts{0, 0, 0, 0};

    // set path prefixes for library files
    io_cfg_type io_cfg{};

    // get instance to primer sequence settings
    primer_cfg_type primer_cfg{};

    // parse options and init io and primer configurators
    options opt(argc, argv, primer_cfg, io_cfg);

    // create FM index if SKIP_IDX not in argument list
    int ret_code;
    if (io_cfg.skip_idx())
    {
        std::cout << "MESSAGE: skip index recomputation" << std::endl;
    }
    else if ((ret_code = fm_index(io_cfg)))
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
    TKLocations locations;

    // directory info needed for genmap's fasta file parser
    TDirectoryInformation directoryInformation;

    // container for fasta header lines
    TSequenceNames sequenceNames;

    // container for fasta sequence lengths
    TSequenceLengths sequenceLengths;

    start = std::chrono::high_resolution_clock::now();
    // compute k-mer mappings
    fm_map(io_cfg, primer_cfg, locations);
    finish = std::chrono::high_resolution_clock::now();
    // duration obj
    //auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(finish - start);
    runtimes[TIMEIT::MAP] += std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();

    // Do not modify or delete STATS lines, since they are captured for statistical analysis
    std::cout << "\nINFO: kmers init = " << locations.size() << std::endl;

    // filter k-mers by frequency and chemical properties
    // TODO: result structure for references and k-mer pairs: candidates/matches
    // vector storing k-mer IDs and their locations, i.e. {TSeq: [(TSeqAccession, TSeqPos)]}
    start = std::chrono::high_resolution_clock::now();
    TReferences references;
    TKmerIDs kmerIDs;
    // TSeqNoMap seqNoMap;
    transform_and_filter(io_cfg, locations, references, kmerIDs, &kmerCounts);
    finish = std::chrono::high_resolution_clock::now();
    runtimes[TIMEIT::FILTER1_TRANSFORM] += std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();

    std::cout << "\nINFO: kmers after filter1 & transform = " << get_num_kmers(kmerIDs) << std::endl;

    // TODO: delete locations
    using TPairList = TPairList<TPair<TCombinePattern<TKmerID, TKmerLength>>>;
    TPairList pairs;

    start = std::chrono::high_resolution_clock::now();

    // template<typename TPairList>
    // void combine(TReferences const & references, TKmerIDs const & kmerIDs, TPairList & pairs, TKmerCounts & stats)

    combine<TPairList>(references, kmerIDs, pairs, &kmerCounts);
    finish = std::chrono::high_resolution_clock::now();
    runtimes[TIMEIT::COMBINE_FILTER2] += std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();

    std::cout << "INFO: pairs combined = " << get_num_pairs<TPairList>(pairs) << std::endl;

    // Decomment the following line for analysing unique kmer combinations.
    start = std::chrono::high_resolution_clock::now();

// template<typename TPairList, typename TPairFreq>
// void filter_pairs(TReferences & references, TKmerIDs const & kmerIDs, TPairList & pairs, std::vector<TPairFreq> & pair_freqs, TKmerCounts * kmerCounts = nullptr)
    // List of unique pair frequencies
    TPairFreqList pair_freqs;

    filter_pairs<TPairList, TPairFreqList>(references, kmerIDs, pairs, pair_freqs, &kmerCounts);

    finish = std::chrono::high_resolution_clock::now();
    runtimes[TIMEIT::PAIR_FREQ] += std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();

    std::cout << "INFO: pairs after frequency cutoff = " << get_num_pairs<TPairList>(pairs) << std::endl;


    //////////////////////////////////////////
    std::cout << "MESSAGE: ... done." << std::endl;

    std::cout << "K\tMAP\t\tFILTER1_TRANSFORM\tCOMBINE_FILTER2\tPAIR_FREQ\t|\tSUM [μs]\n" << std::string(100, '_') << "\n";
    std::cout << "[" << 16 << ":" << 25 << "]\t" << runtimes[priset::TIMEIT::MAP] << "\t" <<
            '\t' << runtimes[priset::TIMEIT::FILTER1_TRANSFORM] <<
            '\t' << runtimes[priset::TIMEIT::COMBINE_FILTER2] << '\t' << runtimes[priset::TIMEIT::PAIR_FREQ] <<
            "\t|\t" << std::accumulate(std::cbegin(runtimes), std::cend(runtimes), 0) << '\n';


    return 0;
}
