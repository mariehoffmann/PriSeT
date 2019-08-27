#include <cassert>
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
*/

/* encode a single sequence (K_min = 0) or a list of sequences with same prefix
 * and lengths in [K_min : seq.size()].
 * bits [0 .. |seq|-1]  complete sequence
 * bit [|seq|]          closure symbol 'C'
 * bits [60:64]         lower sequence length bound in case of variable length
 */
uint64_t dna_encoder_test(priset::TSeq const & seq, TKmerLength const K_min = 0)
{
    uint64_t code = (K_min << 60ULL) + (1ULL << uint64_t(seqan::length(seq) << 1ULL));
    //std::cout << "length K_min = " << K_min << " encoded as: " << code << std::endl;
    for (uint64_t i = 0; i < seqan::length(seq); ++i)
        code += uint64_t(seq[i]) * (1ULL << (i << 1ULL));
    return code;
}

// return full length sequence, ignore variable length info in leading bits.
TSeq dna_decoder_test(uint64_t code)
{
    #define NDEBUG
    assert(code != 0);
    code &= (1ULL << 60ULL) - 1ULL;
    std::array<std::string, 4> sigmas = {"A", "C", "G", "T"};
    priset::TSeq d = "";
    while (code != 1)
    {
        d += sigmas[3 & code];
        code >>= 2;
    }
    return d;
}

void dna_decoder_test(uint64_t code, std::vector<TSeq> & decodes)
{
    // note that assert converted to nop due to seqan's #define NDEBUG
    if (code == 0ULL)
        throw std::invalid_argument("ERROR: invalid argument for decoder, code > 0.");
    decodes.clear();
    // extract 4 highest bits encoding variable kmer length
    uint64_t min_K = ((15ULL << 60ULL) & code) >> 60ULL;
    code &= (1ULL << 60ULL) - 1ULL;
    std::cout << "min_K = " << min_K << std::endl;
    std::array<std::string, 4> sigmas = {"A", "C", "G", "T"};
    priset::TSeq d = "";
    uint64_t k = 1;
    while (code != 1)
    {
        std::cout << "current decode = " << d << ", current code = " << code << std::endl;
        d += sigmas[3 & code];
        if (min_K && min_K <= k++)
            decodes.push_back(d);
        code >>= 2;
    }
}

/*
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

int main()
{
    priset::TSeq seq = "ACGTACGTAAAA";
    std::cout << seq << " as full length code [12] = " << dna_encoder_test(seq) << " and back to " << dna_decoder_test(dna_encoder_test(seq)) << std::endl;
    auto c = dna_encoder_test(seq, 9);
    std::cout << seq << " as variable length code [9:12]= " << c << std::endl;
    std::vector<TSeq> decodes;
    dna_decoder_test(c, decodes);
    std::cout << "decoded: ";
    for (auto decode : decodes)
        std::cout << decode << ", " << std::endl;

    std::cout << "code 0: \n";
    dna_decoder_test(0);
    //converter_test();
}
