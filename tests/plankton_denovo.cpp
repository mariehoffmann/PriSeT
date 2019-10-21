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

// ./clade_X $taxid /Volumes/plastic_data/tactac/subset/$taxid /Volumes/plastic_data/priset/work/$taxid
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
        if (!fs::exists(tmp_dir))
            fs::create_directory(tmp_dir);

    }
};

struct TPrimerKey
{
    uint64_t fwd;
    uint64_t rev;
    TPrimerKey(uint64_t fwd_, uint64_t rev_) : fwd(fwd_), rev(rev_) {}
    bool operator==(TPrimerKey const & lhs) const
    {
        return (fwd == lhs.fwd && rev == lhs.rev);
    }
};

// hash defining in global namespace is ill-formed, needs to be declared in same namespace like std::hash!
//namespace std
//{ // doesn't work with struct hash, overwriting the one in std::hash
struct hash_pp
{

    template<typename TPrimerKey>
    std::size_t operator()(TPrimerKey const & key) const
    {
        return std::hash<std::string>()(std::to_string(key.fwd) + std::to_string(key.rev));
    }
};


int main(int argc, char ** argv)
{
    chemical_filter_test();

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
    fm_map(io_cfg, primer_cfg, locations);
    std::cout << "INFO: kmers init = " << std::accumulate(locations.begin(), locations.end(), 0, [](unsigned ctr, auto & location){return ctr + location.second.first.size();}) << std::endl;

    TReferences references;
    TKmerIDs kmerIDs;
    TSeqNoMap seqNoMap;
    filter_and_transform(io_cfg, locations, references, kmerIDs, seqNoMap, kmerCounts);
    std::cout << "INFO: kmers after filter1 & transform = " << get_num_kmers(kmerIDs) << std::endl;

    // TODO: delete locations
    using TPairList = TPairList<TPair<TCombinePattern<TKmerID, TKmerLength>>>;
    TPairList pairs;
    // dictionary collecting (unique) pair frequencies

    combine<TPairList, TPrimerKey>(references, kmerIDs, pairs, kmerCounts);
    std::cout << "INFO: pairs after combiner = " << get_num_pairs<TPairList>(pairs) << std::endl;

    //  <freq, <kmer_fwd, mask_fwd, kmer_rev, mask_rev> >
    using TPairFreq = std::pair<uint32_t, std::tuple<uint64_t, uint64_t, uint64_t, uint64_t>>;
    std::vector<TPairFreq> pair_freqs;
    filter_pairs(references, kmerIDs, pairs, pair_freqs);

    // Sort pairs by frequency
    std::sort(pair_freqs.begin(), pair_freqs.end(), [&](TPairFreq & p1, TPairFreq & p2){return p1.first > p2.first});
    std::cout << "Frequency\tForward\t\tReverse\n";
    for (auto k = 0; k < std::max(pair_freqs.size(), 10); ++k)
    {
        TPairFreq pf = pair_freqs.at(k);
        std::cout << pf.first << "\t" << dna_decoder(std::get<0>(pf.second), std::get<1>(pf.second)) << "\t";
        std::cout << dna_decoder(std::get<2>(pf.second), std::get<3>(pf.second)) << std::endl; 
    }

    return 0;
}
