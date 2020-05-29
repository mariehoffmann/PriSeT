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

#include "../src/priset.hpp"
#include "../src/types.hpp"

namespace fs = std::experimental::filesystem;

// g++ ../PriSeT/tests/sample_pairs_test.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -lstdc++fs -Wno-unknown-pragmas -lstdc++fs -DNDEBUG -O3 -I/Users/troja/include -L/Users/troja/lib -lsdsl -ldivsufsort -o sample_pairs_test

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

    void tear_down()
    {
        // delete index dir
        if (fs::remove_all(idx_dir))
            std::cout << "ERROR: could not remove idx_dir = " << idx_dir << std::endl;
        // delete tmp dir
        if (fs::remove_all(tmp_dir))
            std::cout << "ERROR: could not remove tmp_dir = " << tmp_dir << std::endl;
    }
};

int main(/*int argc, char ** argv*/)
{
    setup su{};

    unsigned const argc = 6;
    char * const argv[argc] = {"priset", "-l", &su.lib_dir[0], "-w", &su.work_dir[0], "-s"};
    for (unsigned i = 0; i < argc; ++i) std::cout << argv[i] << " ";
    std::cout << std::endl;

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

    TKLocations locations;
    TDirectoryInformation directoryInformation;
    TSequenceNames sequenceNames;
    TSequenceLengths sequenceLengths;

    // compute k-mer mappings
    fm_map(io_cfg, primer_cfg, locations);

    std::cout << "INFO: kmers init = " << locations.size() << std::endl;

    TReferences references;
    TKmerIDs kmerIDs;
    TSeqNoMap seqNoMap;
    transform_and_filter(io_cfg, locations, references, kmerIDs, seqNoMap, kmerCounts);

    std::cout << "INFO: kmers after filter1 & transform = " << get_num_kmers(kmerIDs) << std::endl;

    // TODO: delete locations
    using TPairList = TPairList<TPair<TCombinePattern<TKmerID, TKmerLength>>>;
    TPairList pairs;

    combine<TPairList>(references, kmerIDs, pairs, kmerCounts);

    filter_pairs(references, kmerIDs, pairs);

    std::cout << "INFO: pairs after frequency cutoff = " << get_num_pairs<TPairList>(pairs) << std::endl;

    // sample of 5 randomly chosen pairss
    for (uint8_t i = 0; i < 5; ++i)
    {
        std::srand(std::time(nullptr) + i);
        uint idx = std::rand() % pairs.size();
        std::cout << "sampled pair index: " << idx << std::endl;
        TPair<TCombinePattern<TKmerID, TKmerLength>> pair = pairs.at(idx);
        TKmerID kmer_fwd = kmerIDs[pair.reference][pair.r_fwd];
        TKmerID kmer_rev = kmerIDs[pair.reference][pair.r_rev];
        if (pair.cp.none())
        {
            std::cout << "ERROR: no bit set in combination bitset!\n";
        }
        // sample random combination
        idx = std::rand() % 100;
        while (!pair.cp[idx])
            idx = (idx + 1) % 100;
        std::cout << "Reference ID_cx\tPrimer fwd\tPrimer rev\n";
        std::cout << pair.reference << "\t" << "\t" << dna_decoder(kmer_fwd, ONE_LSHIFT_63 >> (idx/10));
        std::cout << "\t" << dna_decoder(kmer_rev, ONE_LSHIFT_63 >> (idx % 10)) << std::endl;
    }

    su.tear_down();

    return 0;
}
