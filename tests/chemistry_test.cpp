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

#include "../combine_types.hpp"
#include "../filter.hpp"
#include "../primer_cfg_type.hpp"
#include "../types.hpp"
#include "../utilities.hpp"

namespace fs = std::experimental::filesystem;
using namespace priset;

// g++ ../PriSeT/src/tests/chemistry_test.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -lstdc++fs -DNDEBUG -O3 -I/Users/troja/include -L/Users/troja/lib -lsdsl -ldivsufsort -I .. -o chemistry_test

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
    //primer_cfg_type const & primer_cfg{};
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
}

// test Tm and CG content in range, consecutive runs and di-nucleotides are also
// called in corresponding function, but testes by a separate test (see filter_repeats_runs_test).
void chemical_filter_single_pass_test()
{
    // PRIMER_MIN_TM 50.0, PRIMER_MAX_TM 62.0, CG_MIN_CONTENT .4, CG_MAX_CONTENT .6

    // case: all in range, CG = 9 (0.5625), Tm = 7*2 + 9*4 = 50
    TKmerID kmerID = ONE_LSHIFT_63 + dna_encoder("CCGTACGTAACCGGTT");
    TKmerID kmerID_ref = kmerID;
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
}

int main()
{
//    test_dTM();
//    filter_repeats_runs_test();
    filter_CG_clamp_test();
//    chemical_filter_single_pass_test();
    return 0;
}
