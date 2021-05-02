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
#include "../src/algorithm.hpp"
#include "../src/fm.hpp"
#include "../src/IOConfig.hpp"
#include "../src/output.hpp"
#include "../src/PrimerConfig.hpp"
#include "../src/taxonomy.hpp"
#include "../src/types.hpp"
#include "../src/utilities.hpp"

namespace fs = std::experimental::filesystem;
using namespace priset;

// static constexpr bool outputProgress = false;

// -mpopcnt gives speedup of 5, otherwise call __popcountdi2 called for __builtin_popcountll
// g++ ../PriSeT/tests/performance_test.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -lstdc++fs -Wno-unknown-pragmas -lstdc++fs -DNDEBUG -O3 -mpopcnt -I/Users/troja/include -L/Users/troja/lib -lsdsl -ldivsufsort -o performance_test
// ./performance_test $taxid /Volumes/plastic_data/tactac/subset/$taxid /Volumes/plastic_data/priset/work/$taxid
// tst on 304574

struct setup
{
    std::string lib_dir = fs::canonical("../PriSeT/test/library/3041").string();
    std::string work_dir = fs::canonical("../PriSeT/test/work/3041").string();

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
        if (!fs::exists(tmp_dir))
            fs::create_directory(tmp_dir);
    }
};

/* Measure runtime for PriSeT components */
int main(int argc, char ** argv)
{
    if (argc != 4)
    {
        std::cout << "Give taxid, and paths to lib and work dirs.\n";
        exit(-1);
    }
    setup su{argv[2], argv[3]};
    std::string taxid = argv[1];
    unsigned const priset_argc = 6;
    char * const priset_argv[priset_argc] = {"priset", "-l", argv[2], "-w", argv[3], "-s"};
    for (unsigned i = 0; i < priset_argc; ++i) std::cout << priset_argv[i] << " ";
    std::cout << std::endl;

    // Store start and finish times for optional runtime measurements
    std::chrono::time_point<std::chrono::system_clock> start, finish;
    std::array<size_t, TIMEIT::SIZE> runtimes;
    runtimes.fill(0);

    // collect number of kmers or kmer pairs left after relevant processing steps
    TKmerCounts kmerCounts{0, 0, 0, 0};

    // set path prefixes for library files
    IOConfig io_cfg{};

    // get instance to primer sequence settings
    PrimerConfig primer_cfg{};

    // parse options and init io and primer configurators
    options opt(priset_argc, priset_argv, primer_cfg, io_cfg);

    TKLocations locations;
    TDirectoryInformation directoryInformation;
    TSequenceNames sequenceNames;
    TSequenceLengths sequenceLengths;

    // compute k-mer mappings
    start = std::chrono::high_resolution_clock::now();
    fm_map(io_cfg, primer_cfg, locations);
    finish = std::chrono::high_resolution_clock::now();
    runtimes[TIMEIT::MAP] += std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();

    std::cout << "INFO: kmers init = " << std::accumulate(locations.begin(), locations.end(), 0, [](unsigned ctr, auto & location){return ctr + location.second.first.size();}) << std::endl;

    // filter k-mers by frequency and chemical properties
    start = std::chrono::high_resolution_clock::now();
    TReferences references;
    TKmerIDs kmerIDs;
    transform_and_filter(io_cfg, locations, references, kmerIDs, &kmerCounts);
    finish = std::chrono::high_resolution_clock::now();
    runtimes.at(TIMEIT::FILTER1_TRANSFORM) += std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();
    // count k-mers
    for (auto kmerIDs_per_reference : kmerIDs)
    {
        for (TKmerID kmerID : kmerIDs_per_reference)
            kmerCounts[KMER_COUNTS::FILTER1_CNT] += __builtin_popcountll(kmerID >> 54);
    }


    std::cout << "INFO: kmers after filter1 & transform = " << get_num_kmers(kmerIDs) << std::endl;

    using TPairList = TPairList<TPair<CombinePattern>>;
    TPairList pairs;
    // dictionary collecting (unique) pair frequencies

    start = std::chrono::high_resolution_clock::now();
    combine<TPairList>(references, kmerIDs, pairs, &kmerCounts);
    finish = std::chrono::high_resolution_clock::now();
    runtimes.at(TIMEIT::COMBINE_FILTER2) += std::chrono::duration_cast<std::chrono::milliseconds>(finish - start).count();

    std::cout << "INFO: pairs after combiner = " << kmerCounts[KMER_COUNTS::COMBINER_CNT] << std::endl;

    std::vector<TPairFreq> pair_freqs;
    start = std::chrono::high_resolution_clock::now();

    // List of unique pair frequencies
    TPairFreqList pair_freqs;

    filter_pairs<TPairList, TPairFreqList>(references, kmerIDs, pairs, pair_freqs, &kmerCounts);

    finish = std::chrono::high_resolution_clock::now();
    runtimes[TIMEIT::PAIR_FREQ] += std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();

    std::cout << "INFO: pairs after frequency cutoff = " << get_num_pairs<TPairList>(pairs) << std::endl;

    std::cout << "MESSAGE: ... done." << std::endl;

    std::cout << "K\tMAP\t\tFILTER1_TRANSFORM\tCOMBINE_FILTER2\tPAIR_FREQ\t|\tSUM [μs]\n" << std::string(100, '_') << "\n";
    std::cout << "[" << 16 << ":" << 25 << "]\t" << runtimes[priset::TIMEIT::MAP] << "\t" <<
            '\t' << runtimes[priset::TIMEIT::FILTER1_TRANSFORM] <<
            '\t' << runtimes[priset::TIMEIT::COMBINE_FILTER2] << '\t' << runtimes[priset::TIMEIT::PAIR_FREQ] <<
            "\t|\t" << std::accumulate(std::cbegin(runtimes), std::cend(runtimes), 0) << "\n\n";

    // Sort pairs by frequency
    std::sort(pair_freqs.begin(), pair_freqs.end(), [&](TPairFreq & p1, TPairFreq & p2){return p1.first > p2.first;});
    std::stringstream sstream;
    std::cout << "#ID\tForward\tReverse\tFrequency\tTm\tCG\n";

    return 0;
}
