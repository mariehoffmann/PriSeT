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

<<<<<<< HEAD:tests/chemistry_test.cpp
//TODO: this ordering does not work, how to make define directive be found (and that it does not occur too late in translation unit)
// #include "../src/combine_types.hpp"
// #include "../src/filter.hpp"
// #include "../src/primer_cfg_type.hpp"
// #include "../src/types.hpp"
// #include "../src/utilities.hpp"

#include "../src/argument_parser.hpp"
#include "../src/combine_types.hpp"
#include "../src/filter.hpp"
#include "../src/fm.hpp"
#include "../src/io_cfg_type.hpp"
#include "../src/primer_cfg_type.hpp"
||||||| merged common ancestors
#include "../src/combine_types.hpp"
#include "../src/filter.hpp"
#include "../src/primer_cfg_type.hpp"
=======
#include "../src/algorithm.hpp"
#include "../src/PrimerConfig.hpp"
>>>>>>> barcode_miner:test/chemistry_test.cpp
#include "../src/types.hpp"
#include "../src/utilities.hpp"

namespace fs = std::experimental::filesystem;
using namespace priset;

<<<<<<< HEAD:tests/chemistry_test.cpp
// g++ ../PriSeT/tests/chemistry_test.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -lstdc++fs -DNDEBUG -O3 -I/Users/troja/include -L/Users/troja/lib -lsdsl -ldivsufsort -I .. -o chemistry_test

//seqan::String<priset::dna> dna_decoder(uint64_t code, uint64_t const mask);
// called in combiner
void test_dTM()
{
    // 4, 12, 8
    // TKmerID kmerID1 = 1513342761473819536;
    // uint64_t mask1 = ONE_LSHIFT_63 >> 3;
    // TKmerID kmerID2 = 288235407220622179;
    // uint64_t mask2 = ONE_LSHIFT_63 >> 5;
    // float d = dTm(kmerID1, mask1, kmerID2, mask2);
    // if (d != 4)
    // {
    //     std::cout << "ERROR: expect 4, got ";
    //     std::cout << "dTm(" << dna_decoder(kmerID1, mask1) << ", " << dna_decoder(kmerID2, mask2) << ") = " << d << std::endl;
    //     exit(0);
    // }
    // mask1 >>= 2;
    // d = dTm(kmerID1, mask1, kmerID2, mask2);
    // if (d != 12)
    // {
    //     std::cout << "ERROR: expect 12, got ";
    //     std::cout << "dTm(" << dna_decoder(kmerID1, mask1) << ", " << dna_decoder(kmerID2, mask2) << ") = " << d << std::endl;
    //     exit(0);
    // }
    // mask1 >>= 2;
    // d = dTm(kmerID1, mask1, kmerID2, mask2);
    // if (d != 8)
    // {
    //     std::cout << "ERROR: expect 8, got ";
    //     std::cout << "dTm(" << dna_decoder(kmerID1, mask1) << ", " << dna_decoder(kmerID2, mask2) << ") = " << d << std::endl;
    //     exit(0);
    // }
    // std::cout << "SUCCESS\n";

    TKmerID kmerID3 = (1ULL << 61) | dna_encoder("ATTCCAGCTCCAATAGCG");
    // AT = 9, GC = 9
    TKmerID kmerID4 = (1ULL << 59) | dna_encoder("GATTAGATACCATCGTAGTC");
    // AT = 9-12 = -3, CG = 9-8 = 1    -3*2 + 1*4 = -2, abs(-2) = 2
    float d = dTm(kmerID3, 1ULL << 61, kmerID4, 1ULL << 59);
    if (d != 2)
    {
        std::cout << "ERROR: expect 2 Kelvin, got " << d << std::endl;
        exit(0);
    }
    std::cout << "SUCCESS for DIAZ\n";

    // same kmers, but as part of larger ones with various lengths
    // DIAZ fwd contained in kmerID_fwd = 11111100|CTAAACCTCATCATTTAGAGGAAGG
    // DIAZ rev contained in kmerID_fwd = 1111100000|GTAGGTGAACCTGCAGAAGG
    TKmerID kmerID5 = (0b111111ULL << 56) | dna_encoder("ATTCCAGCTCCAATAGCGTTTTTTT");
    TKmerID kmerID6 = (0b11111ULL << 59) | dna_encoder("GATTAGATACCATCGTAGTC");
    d = dTm(kmerID5, 1ULL << 61, kmerID6, 1ULL << 59);
    if (d != 2)
    {
        std::cout << "ERROR: expect 2 Kelvin, got " << d << std::endl;
        exit(0);
    }
    std::cout << "SUCCESS for DIAZ embedded\n";

}

||||||| merged common ancestors
// g++ ../PriSeT/tests/chemistry_test.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -lstdc++fs -DNDEBUG -O3 -I/Users/troja/include -L/Users/troja/lib -lsdsl -ldivsufsort -I .. -o chemistry_test

//seqan::String<priset::dna> dna_decoder(uint64_t code, uint64_t const mask);
// called in combiner
void test_dTM()
{
    // 4, 12, 8
    TKmerID kmerID1 = 1513342761473819536;
    uint64_t mask1 = ONE_LSHIFT_63 >> 3;
    TKmerID kmerID2 = 288235407220622179;
    uint64_t mask2 = ONE_LSHIFT_63 >> 5;
    float d = dTm(kmerID1, mask1, kmerID2, mask2);
    if (d != 4)
    {
        std::cout << "ERROR: expect 4, got ";
        std::cout << "dTm(" << dna_decoder(kmerID1, mask1) << ", " << dna_decoder(kmerID2, mask2) << ") = " << d << std::endl;
        exit(0);
    }
    mask1 >>= 2;
    d = dTm(kmerID1, mask1, kmerID2, mask2);
    if (d != 12)
    {
        std::cout << "ERROR: expect 12, got ";
        std::cout << "dTm(" << dna_decoder(kmerID1, mask1) << ", " << dna_decoder(kmerID2, mask2) << ") = " << d << std::endl;
        exit(0);
    }
    mask1 >>= 2;
    d = dTm(kmerID1, mask1, kmerID2, mask2);
    if (d != 8)
    {
        std::cout << "ERROR: expect 8, got ";
        std::cout << "dTm(" << dna_decoder(kmerID1, mask1) << ", " << dna_decoder(kmerID2, mask2) << ") = " << d << std::endl;
        exit(0);
    }
    std::cout << "SUCCESS\n";
}

=======
// TODO: rewrite as unit test for gtest
>>>>>>> barcode_miner:test/chemistry_test.cpp
void filter_repeats_runs_test()
{
    TKmerID kmerID = ONE_LSHIFT_63 + dna_encoder("ACGTAAAAACGTACGT");
    TKmerID kmerID_ref = kmerID - ONE_LSHIFT_63;
    std::cout << "ACGTAAAAACGTACGT to bits with head bit:\t" << bits2str(kmerID) << std::endl;
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
<<<<<<< HEAD:tests/chemistry_test.cpp
    // chemical_filter_single_pass(kmerID);
    // if (kmerID_ref != kmerID)
    //     std::cout << "ERROR: expect " << kmerID2str(kmerID_ref) << " got " << kmerID2str(kmerID) << std::endl;
    // else
    //     std::cout << "INFO: Success\n";
    // std::cout << "kmerID as bitstr = " << bits2str(kmerID) << std::endl;
    //
    // std::cout << "------------------------------------------\n";
    // // case: only CG content out of range, CG: 10 (0.625), Tm = 2*6 + 4*10 = 62
    // kmerID = ONE_LSHIFT_63 + dna_encoder("CCGACCGTAACCGGTT");
    // kmerID_ref = kmerID - ONE_LSHIFT_63;
    // chemical_filter_single_pass(kmerID);
    // if (kmerID_ref != kmerID)
    //     std::cout << "ERROR: expect " << kmerID2str(kmerID_ref) << " got " << kmerID2str(kmerID) << std::endl;
    // else
    //     std::cout << "INFO: Success\n";
    //
    // std::cout << "------------------------------------------\n";
    // // case: only Tm out of range, CG = 7 (0.4375), Tm = 2*9 + 4*7 = 46
    // kmerID = ONE_LSHIFT_63 + dna_encoder("AACGTAACGTAACGGT");
    // kmerID_ref = kmerID - ONE_LSHIFT_63;
    // chemical_filter_single_pass(kmerID);
    // if (kmerID_ref != kmerID)
    //     std::cout << "ERROR: expect " << kmerID2str(kmerID_ref) << " got " << kmerID2str(kmerID) << std::endl;
    // else
    //     std::cout << "INFO: Success\n";
    //
    // // case: Tm out of range for length 16, not for 17
    // // case: only Tm out of range,
    // // l = 16: CG = 6 (0.375!), Tm = 4*6 + 2*10 = 44!
    // // l = 17: CG = 6 (0.35!), Tm = 4*6 + 2*11 = 46!
    // // l = 25: CG = 8 (0.40), Tm = 2*12 + 4*8 = 56
    // kmerID = ONE_LSHIFT_63 + (ONE_LSHIFT_63 >> 1) + (ONE_LSHIFT_63 >> 4) + dna_encoder("AACGTAACGTACTATCACTC");
    // kmerID_ref = kmerID - ONE_LSHIFT_63 - (ONE_LSHIFT_63 >> 1);
    //
    // chemical_filter_single_pass(kmerID);
    // if (kmerID_ref != kmerID)
    //     std::cout << "ERROR: expect " << kmerID2str(kmerID_ref) << " got " << kmerID2str(kmerID) << std::endl;
    // else
    //     std::cout << "INFO: Success\n";

    //11000000|AGAGGGAGCATGAGAAATG to 0|AGAGGGAGCATGAGAAATG
    // 1100000000|TTATGGTAGAGCTGTAT
    kmerID = (1ULL << 63) | (1ULL << 62) | dna_encoder("TTATGGTAGAGCTGTAT");
    auto kmerID_cp{kmerID};
    chemical_filter_single_pass(kmerID);
    if (kmerID != kmerID_cp)
        std::cout << "ERROR: expect unmodified kmerID, got " << kmerID2str(kmerID) << std::endl;
    else
        std::cout << "INFO: Success\n";


    // kmerID = (1ULL << 63) | (1ULL << 62) | dna_encoder("GGATAGTTGGGGGCATC");
    // chemical_filter_single_pass(kmerID);
    // if (!(kmerID && PREFIX_SELECTOR))
    //     std::cout << "ERROR: expect invalid kmerID, got " << kmerID2str(kmerID) << std::endl;
    // else
    //     std::cout << "INFO: Success\n";

||||||| merged common ancestors
    chemical_filter_single_pass(kmerID);
    if (kmerID_ref != kmerID)
        std::cout << "ERROR: expect " << kmerID2str(kmerID_ref) << " got " << kmerID2str(kmerID) << std::endl;
    else
        std::cout << "INFO: Success\n";
    std::cout << "kmerID as bitstr = " << bits2str(kmerID) << std::endl;

    std::cout << "------------------------------------------\n";
    // case: only CG content out of range, CG: 10 (0.625), Tm = 2*6 + 4*10 = 62
    kmerID = ONE_LSHIFT_63 + dna_encoder("CCGACCGTAACCGGTT");
    kmerID_ref = kmerID - ONE_LSHIFT_63;
    chemical_filter_single_pass(kmerID);
    if (kmerID_ref != kmerID)
        std::cout << "ERROR: expect " << kmerID2str(kmerID_ref) << " got " << kmerID2str(kmerID) << std::endl;
    else
        std::cout << "INFO: Success\n";

    std::cout << "------------------------------------------\n";
    // case: only Tm out of range, CG = 7 (0.4375), Tm = 2*9 + 4*7 = 46
    kmerID = ONE_LSHIFT_63 + dna_encoder("AACGTAACGTAACGGT");
    kmerID_ref = kmerID - ONE_LSHIFT_63;
    chemical_filter_single_pass(kmerID);
    if (kmerID_ref != kmerID)
        std::cout << "ERROR: expect " << kmerID2str(kmerID_ref) << " got " << kmerID2str(kmerID) << std::endl;
    else
        std::cout << "INFO: Success\n";

    // case: Tm out of range for length 16, not for 17
    // case: only Tm out of range,
    // l = 16: CG = 6 (0.375!), Tm = 4*6 + 2*10 = 44!
    // l = 17: CG = 6 (0.35!), Tm = 4*6 + 2*11 = 46!
    // l = 25: CG = 8 (0.40), Tm = 2*12 + 4*8 = 56
    kmerID = ONE_LSHIFT_63 + (ONE_LSHIFT_63 >> 1) + (ONE_LSHIFT_63 >> 4) + dna_encoder("AACGTAACGTACTATCACTC");
    kmerID_ref = kmerID - ONE_LSHIFT_63 - (ONE_LSHIFT_63 >> 1);

    chemical_filter_single_pass(kmerID);
    if (kmerID_ref != kmerID)
        std::cout << "ERROR: expect " << kmerID2str(kmerID_ref) << " got " << kmerID2str(kmerID) << std::endl;
    else
        std::cout << "INFO: Success\n";
=======
    filter_Cs(kmerID);
    if (kmerID_ref != kmerID)
        std::cout << "ERROR: expect " << kmerID2str(kmerID_ref) << " got " << kmerID2str(kmerID) << std::endl;
    else
        std::cout << "INFO: Success\n";
    std::cout << "kmerID as bitstr = " << bits2str(kmerID) << std::endl;
>>>>>>> barcode_miner:test/chemistry_test.cpp
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
    std::cout << "passes WWW filter: " << filter_AT_tail(kmerID, '+', ONE_LSHIFT_63 >> 2) << std::endl;
}

void test3()
{
    TKmerID kmerID = (1ULL << 63) | (1ULL << 62) | (1ULL << 61) | (1ULL << 60) | (1ULL << 59) | (1ULL << 58) | (1ULL << 57) | dna_encoder("TTATGGTAGAGCTGTATATGAA");
    auto kmerID_cp{kmerID};
    chemical_filter_single_pass(kmerID);

    if (kmerID != kmerID_cp)
        std::cout << "OK: expect modified kmerID, got " << kmerID2str(kmerID) << std::endl;
    else
        std::cout << "ERROR\n";
}

int main()
{
<<<<<<< HEAD:tests/chemistry_test.cpp
    filter_CG_clamp_test();
||||||| merged common ancestors
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
//    chemical_filter_single_pass_test();
=======
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
>>>>>>> barcode_miner:test/chemistry_test.cpp
    return 0;
}
