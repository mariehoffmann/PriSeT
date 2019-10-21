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

// Verify that known primers pass the chemical filter
void chemical_filter_test()
{
    std::vector<std::string> primers = {
        "GCGGTAATTCCAGCTCCAATAG",
        "TATTCGTATTCCATTGTCAGAG",
        "CAGCAGCCGCGGTAATTCC",
        "GCTTAATTTGACTCAACACGGG",
        "GCTTGTCTCAAAGATTAAGCC",
        "TCCAAGGAAGGCAGCAGGC",
        "GTACACACCGCCCGTC",
        "GTAGGTGAACCTGCAGAAGGATCA",
        "GGACAAAAAGACCCTATG", "GGATAACAGGCTGATCT", "GGACAGAAAGACCCTATG", "GGATAACAGGCTGATCT",
        "CCAGCACCCGCGGTAATTCC", "CCAGCACCTGCGGTAATTCC", "CCAGCAGCCGCGGTAATTCC", "CCAGCAGCTGCGGTAATTCC",
        "TCAATCAAGAACGAAAGT", "TCGATCAAGAACGAAAGT", "TTAATCAAGAACGAAAGT", "TTGATCAAGAACGAAAGT"};

    for (auto primer : primers)
    {
        uint64_t mask = ONE_LSHIFT_63 >> (primer.size() - PRIMER_MIN_LEN);
        uint64_t code = dna_encoder(primer) + mask;
        chemical_filter_single_pass(code);
        if (PREFIX_SELECTOR & code)
            std::cout << primer << " OK\n";
        else
            std::cout << primer << " DID NOT PASS FILTER!\n";
    }
}

// Note: reverse primers are noted in 5' to 3' direction for the forward sting
// (no reverse complement like needed for real PCR!), because we search on same
// string directly in combiner. So, we have to reverse complement a published primer
// eg. noted as GCTTAATTTGACTCAACACGGG
// For later display the reverse sequence would be transformed.
template<typename TPrimerKey>
void load_primers_known(std::unordered_map<TPrimerKey, std::string, hash_pp> & pairs_known)
{
    pairs_known[TPrimerKey{dna_encoder("GCGGTAATTCCAGCTCCAATAG"), dna_encoder("TATTCGTATTCCATTGTCAGAG")}] = "DIV4";
    pairs_known[TPrimerKey{dna_encoder("CAGCAGCCGCGGTAATTCC"), dna_encoder("GCTTAATTTGACTCAACACGGG")}] = "EUK14";
    pairs_known[TPrimerKey{dna_encoder("GCTTGTCTCAAAGATTAAGCC"), dna_encoder("TCCAAGGAAGGCAGCAGGC")}] = "nSSU";
    // 18411191439941538161, decoded: GTACACACCGCCCGTCGCACCTAC with length bits: 1100
    pairs_known[TPrimerKey{dna_encoder("GTACACACCGCCCGTC"), dna_encoder("GTAGGTGAACCTGCAGAAGGATCA")}] = "V9";
    pairs_known[TPrimerKey{dna_encoder("GGACAAAAAGACCCTATG"), dna_encoder("GGATAACAGGCTGATCT")}] = "23S_1";  // GGACARAAAGACCCTATG
    pairs_known[TPrimerKey{dna_encoder("GGACAGAAAGACCCTATG"), dna_encoder("GGATAACAGGCTGATCT")}] = "23S_2";
    // combine unfolded ones, EUK15_11, i.e. EUK15_fwds[2] - EUK15_revs[2] found
    std::vector<std::string> EUK15_fwds{"CCAGCACCCGCGGTAATTCC", "CCAGCACCTGCGGTAATTCC", "CCAGCAGCCGCGGTAATTCC", "CCAGCAGCTGCGGTAATTCC"};
    std::vector<std::string> EUK15_revs{"TCAATCAAGAACGAAAGT", "TCGATCAAGAACGAAAGT", "TTAATCAAGAACGAAAGT", "TTGATCAAGAACGAAAGT"}; //TYRATCAAGAACGAAAGT

    uint16_t i = 1;
    for (std::string EUK15_fwd : EUK15_fwds)
    {
        for (std::string EUK15_rev : EUK15_revs)
            pairs_known[TPrimerKey{dna_encoder(EUK15_fwd), dna_encoder(EUK15_rev)}] = "EUK15_" + std::to_string(i++);
    }
}

// same as combine, but accepts primer pair set to be verified
template<typename TPairList, typename TPrimerKey>
void combine2(TReferences const & references, TKmerIDs const & kmerIDs, TPairList & pairs, TKmerCounts & stats, std::unordered_map<uint64_t, uint32_t> & pair2freq, std::unordered_map<TPrimerKey, std::string, hash_pp> & pairs_known, std::unordered_set<std::string> & verified)
{
    pairs.clear();
    auto cp_ctr = 0ULL;
    for (uint64_t i = 0; i < references.size(); ++i)
    {
        sdsl::bit_vector reference;
        sdsl::util::assign(reference, references.at(i));
        sdsl::rank_support_v5<1, 1> r1s(&references.at(i)); // replace after bugfix with
        sdsl::select_support_mcl<1> s1s(&reference);

        for (uint64_t r_fwd = 1; r_fwd < r1s.rank(reference.size()); ++r_fwd)
        {
            uint64_t idx_fwd = s1s.select(r_fwd);  // text position of r-th k-mer
            TKmerID const kmerID_fwd = kmerIDs[i][r_fwd - 1];
            if (!(kmerID_fwd >> CODE_SIZE))
            {
                std::cout << "ERROR: k length pattern is zero\n";
                exit(-1);
            }
            uint64_t w_begin = idx_fwd + PRIMER_MIN_LEN + TRANSCRIPT_MIN_LEN;
            uint64_t w_end = std::min(reference.size(), idx_fwd + PRIMER_MAX_LEN + TRANSCRIPT_MAX_LEN + 1); // + 1: rank excludes upper bound
            for (uint64_t r_rev = r1s.rank(w_begin) + 1; r_rev <= r1s.rank(w_end); ++r_rev)
            {
                TCombinePattern<TKmerID, TKmerLength> cp;
                uint64_t mask_fwd = ONE_LSHIFT_63;
                while ((((mask_fwd - 1) << 1) & kmerID_fwd) >> 54)
                {
                    if ((mask_fwd & kmerID_fwd)) // && filter_CG_clamp(kmerID_fwd, '+', mask_fwd))
                    {
                        TKmerID const kmerID_rev = kmerIDs.at(i).at(r_rev - 1);
                        uint64_t mask_rev = ONE_LSHIFT_63;

                        while ((((mask_rev - 1) << 1) & kmerID_rev) >> 54)
                        {
                            bool flag = false;
                            cp_ctr++;
                            if (cp_ctr == 1ULL << 63)
                            {
                                std::cout << "cp_ctr reached 1 << 63, exit\n";
                                exit(0);
                            }
                            if ((mask_rev & kmerID_rev) ) //&& filter_CG_clamp(kmerID_rev, '-'))
                            {
                                if (dTm(kmerID_fwd, mask_fwd, kmerID_rev, mask_rev) <= PRIMER_DTM)
                                {
                                    auto code_fwd = get_code(kmerID_fwd, mask_fwd);
                                    auto code_rev = get_code(kmerID_rev, mask_rev);
                                    auto h = hash_pair(code_fwd, code_rev);
                                    if (pair2freq.find(h) == pair2freq.end())
                                        pair2freq[h] = 1;
                                    else
                                        pair2freq[h]++;
                                    cp.set(mask_fwd, mask_rev);
                                    TPrimerKey key{get_code(kmerID_fwd, mask_fwd), get_code(kmerID_rev, mask_rev)};
                                    if (pairs_known.find(key) != pairs_known.end())
                                    {
                                        std::cout << "INFO: primer pair <" << pairs_known[key] << "> found for refID = " << i << " at " << s1s.select(r_fwd) << " and " << s1s.select(r_rev) << std::endl;
                                        verified.insert(pairs_known[key]);
                                    }
                                    ++stats[KMER_COUNTS::COMBINER_CNT];
                                    flag = false;
                                }
                                flag = false;
                            }
                            mask_rev >>= 1; // does not affect search window, since starting position is fixed
                        } // length mask_rev
                    }
                    mask_fwd >>= 1;
                } // length mask_fwd

                if (cp.is_set())
                    pairs.push_back(TPair<TCombinePattern<TKmerID, TKmerLength>>{i, r_fwd, r_rev, cp});
            } // kmerID rev
        } // kmerID fwd
    }
}


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
    std::unordered_map<uint64_t, uint32_t> pair2freq;
    std::unordered_map<TPrimerKey, std::string, hash_pp> pairs_known;
    load_primers_known<TPrimerKey>(pairs_known);
    std::unordered_set<std::string> verified;

    combine2<TPairList, TPrimerKey>(references, kmerIDs, pairs, kmerCounts, pair2freq, pairs_known, verified);
    std::cout << "INFO: pairs after combiner = " << get_num_pairs<TPairList>(pairs) << std::endl;

    std::cout << "Verified primers for current clade: \n";
    if (!verified.size())
        std::cout << 0 << std::endl;
    else
    {
        for (std::string pp : verified)
            std::cout << pp << std::endl;
    }

    return 0;
}
