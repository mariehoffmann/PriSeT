#include <algorithm>
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

// g++ ../PriSeT/tests/plankton_denovo.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -lstdc++fs -Wno-unknown-pragmas -lstdc++fs -DNDEBUG -O3 -I/Users/troja/include -L/Users/troja/lib -lsdsl -ldivsufsort -o denovo

// ./denovo $taxid /Volumes/plastic_data/tactac/subset/$taxid /Volumes/plastic_data/priset/work/$taxid

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

    // timing
    std::chrono::time_point<std::chrono::system_clock> start, finish;
    std::array<size_t, TIMEIT::SIZE> runtimes;

    // collect number of kmers or kmer pairs left after relevant processing steps
    TKmerCounts kmerCounts{0, 0, 0, 0};

    // set path prefixes for library files
    io_cfg_type io_cfg{};

    // get instance to primer sequence settings
    primer_cfg_type primer_cfg{};

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
    runtimes.at(TIMEIT::MAP) += std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();

    std::cout << "INFO: kmers init = " << std::accumulate(locations.begin(), locations.end(), 0, [](unsigned ctr, auto & location){return ctr + location.second.first.size();}) << std::endl;

    TReferences references;
    TKmerIDs kmerIDs;
    TSeqNoMap seqNoMap;
    start = std::chrono::high_resolution_clock::now();
    filter_and_transform(io_cfg, locations, references, kmerIDs, seqNoMap, kmerCounts);
    finish = std::chrono::high_resolution_clock::now();
    runtimes.at(TIMEIT::FILTER1_TRANSFORM) += std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();
    std::cout << "INFO: kmers after filter1 & transform = " << get_num_kmers(kmerIDs) << std::endl;

    // TODO: delete locations
    using TPairList = TPairList<TPair<TCombinePattern<TKmerID, TKmerLength>>>;
    TPairList pairs;
    // dictionary collecting (unique) pair frequencies

    //void combine(TReferences const & references, TKmerIDs const & kmerIDs, TPairList & pairs, TKmerCounts & stats)
    start = std::chrono::high_resolution_clock::now();
    combine(references, kmerIDs, pairs, kmerCounts);
    finish = std::chrono::high_resolution_clock::now();
    runtimes.at(TIMEIT::COMBINE_FILTER2) += std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();

    std::cout << "INFO: pairs after combiner = " << kmerCounts[KMER_COUNTS::COMBINER_CNT] << std::endl;

    //  <freq, <kmer_fwd, mask_fwd, kmer_rev, mask_rev> >
    using TPairFreq = std::pair<uint32_t, std::tuple<uint64_t, uint64_t, uint64_t, uint64_t>>;
    std::vector<TPairFreq> pair_freqs;
    start = std::chrono::high_resolution_clock::now();
    filter_pairs(references, kmerIDs, pairs, pair_freqs, kmerCounts);
    std::cout << "INFO: pairs after pair_freq filter = " << kmerCounts[KMER_COUNTS::FILTER2_CNT] << std::endl;
    finish = std::chrono::high_resolution_clock::now();
    runtimes.at(TIMEIT::PAIR_FREQ) += std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count();

    std::cout << "K\tMAP\t\tFILTER1_TRANSFORM\tCOMBINE_FILTER2\tPAIR_FREQ\t|\tSUM [μs]\n" << std::string(100, '_') << "\n";
    std::cout << "[" << 16 << ":" << 25 << "]\t" << runtimes[priset::TIMEIT::MAP] << "\t" <<
            '\t' << runtimes[priset::TIMEIT::FILTER1_TRANSFORM] <<
            '\t' << runtimes[priset::TIMEIT::COMBINE_FILTER2] << '\t' << runtimes[priset::TIMEIT::PAIR_FREQ] <<
            "\t|\t" << std::accumulate(std::cbegin(runtimes), std::cend(runtimes), 0) << '\n';

    // Sort pairs by frequency
    std::sort(pair_freqs.begin(), pair_freqs.end(), [&](TPairFreq & p1, TPairFreq & p2){return p1.first > p2.first;});
    std::cout << "ID\tFrequency\tForward\tReverse\tTm\tCG\n";
    for (size_t k = 0; k < std::min(pair_freqs.size(), size_t(15)); ++k)
    {
        TPairFreq pf = pair_freqs.at(k);
        const auto & [code_fwd, mask_fwd, code_rev, mask_rev] = pf.second;
        std::string fwd = dna_decoder(code_fwd, mask_fwd);
        std::string rev = reverse_complement(dna_decoder(code_rev, mask_rev));
        std::stringstream sstream;
        sstream << std::hex << std::hash<std::string>()(fwd + rev);
        std::cout << sstream.str().substr(0,8) << "," << pf.first << "," << fwd << ",";
        std::cout << rev << ",\"[" << float(Tm(code_fwd, mask_fwd)) << ",";
        std::cout << float(Tm(code_rev, mask_rev)) << "]\",\"[" << CG(code_fwd, mask_fwd) << "," << CG(code_rev, mask_rev) << "]\"" << std::endl;
    }

    return 0;
}
