#include <cassert>
#include <chrono>
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

#include "../src/filter.hpp"
#include "../src/PrimerConfig.hpp"
#include "../src/types.hpp"
#include "../src/utilities.hpp"

namespace fs = std::experimental::filesystem;
using namespace priset;


void filter_repeats_runs_test()
{
    // |ACGTAAAAACGTACGT| = 16   0000000000000000000000000000000100011011000000000001101100011011

    TKmerID kmerID = ONE_LSHIFT_63 + dna_encoder("ACGTAAAAACGTACGT");
    TKmerID kmerID_ref = kmerID - ONE_LSHIFT_63;
    std::cout << "ACGTAAAAACGTACGT to bits with head bit:\t" << bits2str(kmerID) << std::endl;
    // head should be 0 after return
    filter_repeats_runs(kmerID);
    std::cout << "ACGTAAAAACGTACGT to bits after head rem:\t" << bits2str(kmerID) << std::endl;

    if (kmerID_ref != kmerID)
    {
        std::cout << "ERROR: expect 0 header, got " << bits2str((PREFIX_SELECTOR & kmerID) >> 54) << std::endl;
        exit(0);
    }
    // |TTTTTCGTAAAAGACGTACGT| = 21
    kmerID = (ONE_LSHIFT_63 >> 5) + dna_encoder("TTTTTCGTAAAAGACGTACGT");
    kmerID_ref = kmerID - (ONE_LSHIFT_63 >> 5);
    filter_repeats_runs(kmerID);
    if (kmerID_ref != kmerID)
    {
        std::cout << "ERROR: expect zero header, got: " << bits2str((PREFIX_SELECTOR & kmerID) >> 54) << std::endl;
        exit(0);
    }
    std::cout << "INFO: Success, 5Xs filtered\n";

    // case: multiple lengths, only largest filtered out |ACGTACGTACGTAAAA| = 16, |ACGTACGTACGTAAAAATTTT| = 21
    kmerID = ONE_LSHIFT_63 + (ONE_LSHIFT_63 >> 5) + dna_encoder("ACGTACGTACGTAAAAATTTT");
    kmerID_ref = ONE_LSHIFT_63 + dna_encoder("ACGTACGTACGTAAAAATTTT");
    filter_repeats_runs(kmerID);
    if (kmerID != kmerID_ref)
    {
        std::cout << "ERROR: kmerID_f is 1000000000 ? " << bits2str((PREFIX_SELECTOR & kmerID) >> 54);
        std::cout << " or tail different\n";
    }
    else
        std::cout << "INFO: Success, 5As filtered\n";

    // case: multiple lengths, only largest filtered out |ACGTACGTACGTATATATAT| = 20, |ACGTACGTACGTATATATATA| = 21
    kmerID = (ONE_LSHIFT_63 >> 4) + (ONE_LSHIFT_63 >> 5) + dna_encoder("ACGTACGTACGTATATATATA");
    kmerID_ref = kmerID - (ONE_LSHIFT_63 >> 5);
    //std::cout << "kmerID in: " << std::bitset<64>(kmerID) << " as string " << kmerID2str(kmerID) << std::endl;
    filter_repeats_runs(kmerID);

    if (kmerID_ref != kmerID)
        std::cout << "ERROR: kmerID_f head should be 0b0000100000, but is: " << bits2str((PREFIX_SELECTOR & kmerID) >> 54) << std::endl;
    else
        std::cout << "INFO: Success, 5TAs filtered\n";

    kmerID = ONE_LSHIFT_63 + (ONE_LSHIFT_63 >> 7) + dna_encoder("CGCGCGCGCGACGTACGTACGTACGT"); // 26
    kmerID_ref = (kmerID - ONE_LSHIFT_63 - (ONE_LSHIFT_63 >> 7)) >> 6; // result trimmed to initially encoded length
    filter_repeats_runs(kmerID);
    if (kmerID_ref != kmerID)
        std::cout << "ERROR: kmerID_f head should be 0b0000100000, but is: " << bits2str((PREFIX_SELECTOR & kmerID) >> 54) << std::endl;
    else
        std::cout << "INFO: Success, 5TAs filtered\n";
}

void filter_CG_clamp_test()
{
    //PrimerConfig const & primer_cfg{};
    TKmerID kmerID = ONE_LSHIFT_63 + (ONE_LSHIFT_63 >> 1) + dna_encoder("CACGTACGTAACCGGTT"); // l = 17
    // '+' full length CG_ctr = 3
    bool res = filter_CG_clamp(kmerID, '+', ONE_LSHIFT_63 >> 1);
    if (!res)
        std::cout << "ERROR: expect true, got false.\n";
    else
        std::cout << "INFO: Success\n";
    // '+' length 16 =>  CG_ctr = 4
    res = filter_CG_clamp(kmerID, '+', ONE_LSHIFT_63);
    if (res)
        std::cout << "ERROR: expect false, got true.\n";
    else
        std::cout << "INFO: Success\n";
    // '-', CG_ctr = 3
    kmerID = ONE_LSHIFT_63 + (ONE_LSHIFT_63 >> 5) + dna_encoder("CACGTACGTAACCGGTTACGT");
    res = filter_CG_clamp(kmerID, '-');
    if (!res)
        std::cout << "ERROR: expect true, got false.\n";
    else
        std::cout << "INFO: Success\n";
    kmerID = ONE_LSHIFT_63 + (ONE_LSHIFT_63 >> 5) + dna_encoder("CACGCACGTAACCGGTTACGT");
    res = filter_CG_clamp(kmerID, '-');
    if (res)
        std::cout << "ERROR: expect false, got true.\n";
    else
        std::cout << "INFO: Success\n";
    std::cout << "---------------------------------------------\n";
    kmerID = dna_encoder("TTTGACTCAACACGGGGA") + (ONE_LSHIFT_63 >> 1) + (ONE_LSHIFT_63 >> 2);
    res = filter_CG_clamp(kmerID, '+', ONE_LSHIFT_63 >> 1);
    if (res)
        std::cout << "ERROR: expect false, got true.\n";
    else
        std::cout << "INFO: Success\n";

    // now test kmer[0:16] which should pass test
    res = filter_CG_clamp(kmerID, '+', ONE_LSHIFT_63 >> 2);
    if (res)
        std::cout << "ERROR: expect false, got true.\n";
    else
        std::cout << "INFO: Success\n";
}

// test Tm and CG content in range, consecutive runs and di-nucleotides are also
// called in corresponding function, but testes by a separate test (see filter_repeats_runs_test).
void filter_Cs_test()
{
    // PRIMER_MIN_TM 50.0, PRIMER_MAX_TM 62.0, CG_MIN_CONTENT .4, CG_MAX_CONTENT .6

    // case: all in range, CG = 9 (0.5625), Tm = 7*2 + 9*4 = 50
    TKmerID kmerID = ONE_LSHIFT_63 + dna_encoder("CCGTACGTAACCGGTT");
    TKmerID kmerID_ref = kmerID;
    filter_Cs(kmerID);
    if (kmerID_ref != kmerID)
        std::cout << "ERROR: expect " << kmerID2str(kmerID_ref) << " got " << kmerID2str(kmerID) << std::endl;
    else
        std::cout << "INFO: Success\n";
    std::cout << "kmerID as bitstr = " << bits2str(kmerID) << std::endl;

    std::cout << "------------------------------------------\n";
    // case: only CG content out of range, CG: 10 (0.625), Tm = 2*6 + 4*10 = 62
    kmerID = ONE_LSHIFT_63 + dna_encoder("CCGACCGTAACCGGTT");
    kmerID_ref = kmerID - ONE_LSHIFT_63;
    filter_Cs(kmerID);
    if (kmerID_ref != kmerID)
        std::cout << "ERROR: expect " << kmerID2str(kmerID_ref) << " got " << kmerID2str(kmerID) << std::endl;
    else
        std::cout << "INFO: Success\n";
}

void chemical_debug()
{
    TKmerID kmerID1 = 18411191439941538145ULL;
    //, decoded: GTACACACCGCCCGTCGCACCGAC with length bits: 1111111110
    TKmerID V9_fwd = dna_encoder("GTACACACCGCCCGTCGCACCGAC") + PREFIX_SELECTOR - (1ULL << 54);
    std::cout << "reconstructed: " << (V9_fwd == kmerID1) << std::endl;
    filter_Cs(V9_fwd);
    std::cout << "V9 fwd after chem filter: " << kmerID2str(V9_fwd) << std::endl;
    TKmerID V9_rev = dna_encoder("GTAGGTGAACCTGCAGAAGGATCA") +  (1ULL << 55);
    filter_Cs(V9_rev);
    std::cout << "V9 rev after chem filter: " << kmerID2str(V9_rev) << std::endl;

    std::unordered_set<uint64_t> EUK15_fwds{dna_encoder("CCAGCACCCGCGGTAATTCC"), dna_encoder("CCAGCACCTGCGGTAATTCC"), dna_encoder("CCAGCAGCCGCGGTAATTCC"), dna_encoder("CCAGCAGCTGCGGTAATTCC")};
    for (uint64_t code : EUK15_fwds)
    {
        code += 1ULL << 59;
        filter_Cs(code);
        std::cout << "EUK15 fwd variant passed chemical filter: " << ((code & (1ULL << 59)) ? 1 : 0) << std::endl;
    }
    std::unordered_set<uint64_t> EUK15_revs{dna_encoder("TCAATCAAGAACGAAAGT"), dna_encoder("TCGATCAAGAACGAAAGT"), dna_encoder("TTAATCAAGAACGAAAGT"), dna_encoder("TTGATCAAGAACGAAAGT")}; //
    for (uint64_t code : EUK15_revs)
    {
        code += 1ULL << 61;
        filter_Cs(code);
        std::cout << "EUK15 rev variant passed chemical filter: " << ((code & (1ULL << 61)) ? 1 : 0) << std::endl;
    }

    //
    uint64_t mask_rev = 1ULL << 61;
    uint64_t kmerID_rev = dna_encoder("TTAATCAAGAACGAAAGT") + (1ULL << 61); // 46

    uint64_t mask_fwd = 1ULL << 59;
    uint64_t kmerID_fwd = dna_encoder("CCAGCAGCCGCGGTAATTCC") + (1ULL << 59); // 66
    std::cout << "Would pass dTM: ";
    std::cout << (dTm(kmerID_fwd, mask_fwd, kmerID_rev, mask_rev) <= PRIMER_DTM) << std::endl;
    std::cout << "dTm = " << dTm(kmerID_fwd, mask_fwd, kmerID_rev, mask_rev)  << " and PRIMER_DTM = " << PRIMER_DTM << std::endl;
}

void test_filter_WWW()
{
    TKmerID kmerID = 4036873631493718478;
    std::cout << kmerID2str(kmerID) << std::endl;
    std::cout << "passes WWW filter: " << filter_WWW_tail(kmerID, '+', ONE_LSHIFT_63 >> 2) << std::endl;
}

int main()
{
//    reverse_complement_test();
//    complement_test();
    test_filter_self_annealing();
//    test_filter_cross_annealing();
//    test_filter_WWW();
//    test_reverse_complement();
//    chemical_debug();
//    test_dTM();
//    filter_repeats_runs_test();
//    filter_CG_clamp_test();
//    filter_Cs_test();
    return 0;
}
