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

// g++ ../PriSeT/src/test.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -lstdc++fs -DNDEBUG -O3 -I/Users/troja/include -L/Users/troja/lib -lsdsl -ldivsufsort -o test

struct setup
{
    // TODO: make this runnable with arbitrarily located build folders
    std::string lib_dir = (fs::canonical("../PriSeT/src/tests/library/3041")).string();
    std::string work_dir = (fs::canonical("../PriSeT/src/tests/work/3041")).string();
    priset::io_cfg_type io_cfg{};
    priset::primer_cfg_type primer_cfg{};
    priset::TKmerLocations kmer_locations;

    setup()
    {
        // basic init
        int argc = 5;
        char * argv[5] = {"priset", "-l", &lib_dir[0], "-w", &work_dir[0]};
        for (auto i = 0; i < argc; ++i)
            std::cout << "arg[" << i << "] = " << argv[i] << std::endl;

        priset::options opt(argc, argv, primer_cfg, io_cfg);

        // change primer and transcript length ranges
        // k1: [(1,2), (1,75)], i.e. kmer1 occurs in  sequence 1 at position 2 and 75
        // k2: [(1,5), (1,80)], i.e. kmer2 occurs in sequence 2 at positions 20 and 80
        // -k1[2]-k2[20]-------------k1[75]/k2[80]
        primer_cfg.set_primer_length_range(4, 8);
        primer_cfg.set_transcript_range(50, 800);

        priset::TLocation loc1_kmer1{1, 2};
        priset::TLocation loc2_kmer1{1, 75};
        priset::TLocation loc1_kmer2{1, 5};
        priset::TLocation loc2_kmer2{1, 80};

        // init locations
        priset::TKmerLocation::TLocationVec loc_vec1{{loc1_kmer1, loc2_kmer1}};
        priset::TKmerLocation::TLocationVec loc_vec2{{loc1_kmer2, loc2_kmer2}};
        //kmer_locations.push_back(priset::TKmerLocation(1, loc_vec1));
        //kmer_locations.push_back(priset::TKmerLocation(2, loc_vec2));

    }
};

/*
void gui_test()
{
    setup up{};
    if (priset::gui::generate_app(up.io_cfg) && priset::gui::compile_app(up.io_cfg))
        std::cout << "success" << std::endl;
    else
        std::cout << "failed to generate and compile app" << std::endl;
}

void combine_test()
{
    setup s{};

    priset::TKmerPairs pairs{};
    priset::combine(s.primer_cfg, s.kmer_locations, pairs);
    priset::print_pairs(pairs);
}

void lookup_sequences_test()
{
    seqan::StringSet<seqan::DnaString, seqan::Owner<seqan::ConcatDirect<>>> text;
    fs::path text_path = "/Users/troja/priset/335928/work/index/index.txt"; //io_cfg.get_index_txt_path();
    std::cout << "text_path: " << text_path << std::endl;
    seqan::open(text, text_path.string().c_str(), seqan::OPEN_RDONLY);

    std::cout << seqan::valueById(text, 0) << std::endl;
    typedef seqan::Iterator<seqan::StringSet<seqan::DnaString, seqan::Owner<seqan::ConcatDirect<>>>>::Type TStringSetIterator;
    for (TStringSetIterator it = seqan::begin(text); it != seqan::end(text); ++it)
        std::cout << *it << '\n';
}

void dimerization_test()
{

}

uint64_t dna_encoder(priset::TSeq const & seq)
{
    uint64_t code = 0;
    for (uint64_t i = 0; i < seqan::length(seq); ++i)
    {
        code += uint64_t(seq[i]) * (1 << (i << 1));
    }
    std::cout << "sum = " << code << std::endl;
    return code + (uint64_t(1) << uint64_t(seqan::length(seq) << 1));
}

priset::TSeq dna_decoder(uint64_t code)
{
    assert(code > 0);
    std::array<std::string, 4> decodes = {"A", "C", "G", "T"};
    priset::TSeq d = "";
    while (code != 1)
    {

        d += decodes[3 & code];
        code >>= 2;
    }
    return d;
}

void filter_repeats_runs_test()
{
    priset::TSeq seqseq = "GATATATATG";
    uint64_t seq = dna_encoder(seqseq);
    std::cout << "kmer seq = " << seqseq << " passes repeat_and_runs filter: " << priset::filter_repeats_runs(seq) << std::endl;
    seqseq = "GATATATAGATGG";
    seq = dna_encoder(seqseq);
    std::cout << "kmer seq = " << seqseq << " passes repeat_and_runs filter: " << priset::filter_repeats_runs(seq) << std::endl;
    seqseq = "GGGATATATAT";
    seq = dna_encoder(seqseq);
    std::cout << "kmer seq = " << seqseq << " passes repeat_and_runs filter: " << priset::filter_repeats_runs(seq) << std::endl;
    seqseq = "GATATATTTTTGG";
    seq = dna_encoder(seqseq);
    std::cout << "kmer seq = " << seqseq << " passes repeat_and_runs filter: " << priset::filter_repeats_runs(seq) << std::endl;
    seqseq = "GATATAAAAA";
    seq = dna_encoder(seqseq);
    std::cout << "kmer seq = " << seqseq << " passes repeat_and_runs filter: " << priset::filter_repeats_runs(seq) << std::endl;
}
*/

/*
 * seq 0:   00100100000000010010000000000000000000001
 *          01234567890123456789012345678901234567890
 *          kmer1, kmer2, kmer1, kmer3, kmer4
 * seq 2:   01100000000000000100000000000000000000001
 *          01234567890123456789012345678901234567890
 *          kmer1, kmer4, kmer3
 * kmer1 =
*/
void converter_test()
{
    setup su{};
    std::vector<priset::TSeq> kmers = {
            "CACGATTACCAATCAC",
            "GATTACCAATCACGAT",
            "GATTACCAATCACGGC",
            "ACGATTACCAATCACG"};
    for (auto const kmer : kmers)
        std::cout << kmer << " => " << priset::dna_encoder(kmer) << std::endl;

    std::vector<priset::TKmerID> kmer_IDs;
    std::transform(kmers.begin(), kmers.end(), std::back_inserter(kmer_IDs), [](priset::TSeq seq){return priset::dna_encoder(seq);});

    // src
    priset::TKLocations locations;
    std::vector<priset::TLocation> dummy;
    std::vector<priset::TKLocation> keys{priset::TKLocation{0, 2, 16}, priset::TKLocation{0, 5, 16}, priset::TKLocation{0, 18, 16}, TKLocation{0, 50, 16}};
    std::vector<std::vector<TLocation>> values{
        std::vector<TLocation>{TLocation{0,2}, TLocation{0, 15}, TLocation{2, 1}},
        std::vector<TLocation>{TLocation{0,5}, TLocation{2, 50}},
        std::vector<TLocation>{TLocation{0,18}, TLocation{2, 17}},
        std::vector<TLocation>{TLocation{0,50}, TLocation{2, 2}}
    };

    for (unsigned i = 0; i < keys.size(); ++i)
    {
        auto value = std::make_pair(values[i], dummy);
        std::cout << "value.first.size = " << value.first.size() << std::endl;
        locations.insert({keys[i], value});
    }
    // dst
    priset::TReferences references;
    priset::TKmerIDs kmerIDs;
    priset::TSeqNoMap seqNoMap;
    priset::TSeqNo const cutoff = 2;

    std::unordered_map<std::string, priset::TSeq> seq_map;
    seq_map.insert({"0_2",  "CACGATTACCAATCAC"});
    seq_map.insert({"0_5",  "GATTACCAATCACGAT"});
    seq_map.insert({"0_18", "GATTACCAATCACGGC"});
    seq_map.insert({"0_50", "ACGATTACCAATCACG"});

    priset::frequency_filter2(su.io_cfg, locations, references, kmerIDs, seqNoMap, cutoff, seq_map);
    std::cout << "references: \n";
    for (auto reference : references)
    {
        for (unsigned i = 0; i < reference.size(); ++i)
            std::cout << reference[i];
        std::cout << std::endl;
    }
    std::cout << "and associated kmers: \n";
    for (unsigned i = 0; i < kmerIDs.size(); ++i)
    {
        for (priset::TKmerID kmerID : kmerIDs[i])
            std::cout << kmerID << ", ";
        std::cout << std::endl;
    }

}

int main()
{
    converter_test();
}
