#include <array>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <experimental/filesystem>
#include <fstream>
#include <numeric>
#include <regex>
#include <sstream>
#include <string>
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

// g++ ../PriSeT/tests/clade_X.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -lstdc++fs -Wno-unknown-pragmas -lstdc++fs -DNDEBUG -O3 -I/Users/troja/include -L/Users/troja/lib -lsdsl -ldivsufsort -o clade_X
// ./clade_X <taxid> /Volumes/plastic_data/tactac/subset/<taxid> /Volumes/plastic_data/priset/work/<taxid>
struct setup
{
    std::string lib_dir;
    std::string work_dir;

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

template<typename TMatches>
void init_primer_search(fs::path & primer_file, TMatches & matches)
{
    std::ifstream infile(primer_file);
    std::string line;
    std::pair<std::string, std::string> key;
    while (std::getline(infile, line))
    {
        if (line.compare(0, 1, ">"))
        {
            std::cout << "ERROR: expect line starting with '>'\n";
            break;
        }
        auto header = line;
        if (!std::getline(infile, line))
        {
            std::cout << "ERROR: unexpected end of file\n";
            break;
        }
        key = std::pair<std::string, std::string>{header.substr(1, header.size() - 1), line};
        matches[key] = 0;
    }
    std::cout << "matches.size() = " << matches.size() << std::endl;
}

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
    std::cout << "Runtime fm_map: " << std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count() << std::endl;

    std::cout << "INFO: kmers init = " << std::accumulate(locations.begin(), locations.end(), 0, [](unsigned ctr, auto & location){return ctr + location.second.first.size();}) << std::endl;

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
    // dictionary collecting (unique) pair frequencies
    std::unordered_map<uint64_t, uint32_t> pair2freq;

    start = std::chrono::high_resolution_clock::now();
    combine<TPairList>(references, kmerIDs, pairs, kmerCounts, pair2freq);
    finish = std::chrono::high_resolution_clock::now();
    std::cout << "Runtime combine: " << std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count()  << std::endl;

    std::cout << "INFO: pairs after combiner = " << get_num_pairs<TPairList>(pairs) << std::endl;

    start = std::chrono::high_resolution_clock::now();

    filter_pairs(references, kmerIDs, pairs, pair2freq);

    finish = std::chrono::high_resolution_clock::now();
    std::cout << "Runtime filter_pairs: " << std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count() << std::endl;

    std::cout << "INFO: pairs after frequency cutoff = " << get_num_pairs<TPairList>(pairs) << std::endl;

    // get unique primers that are part of pairs into file
    fs::path primer_file = su.tmp_dir / ("primers_" + taxid + ".csv");
    std::unordered_set<std::string> kmers_unique;
    write_primer_file<TKmerIDs, TPairList>(kmerIDs, pairs, primer_file, kmers_unique);

    // file containing published primers as regular expression due to some ambiguous bps
    fs::path primers_known = "/Users/troja/git/PriSet_git2/PriSeT/tests/ecology_primers_regex.fasta";
    std::map<std::pair<std::string, std::string>, uint16_t> matches;
    init_primer_search(primers_known, matches);

    // Iterate over regex and unique kmer strings and report matches
    for (auto it = matches.begin(); it != matches.end(); ++it)
    {
        for (std::string const & kmer_str : kmers_unique)
        {
            std::regex const primer_rx(it->first.second);
            bool found = std::regex_match(kmer_str, primer_rx);
            if (found)
                ++(it->second);
        }
        std::cout << "(" << it->first.first << ", " << it->first.second << "): " << it->second << std::endl;
    }




    return 0;
}
