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

void test_filter_repeats_runs()
{
    // |ACGTAAAAACGTACGT| = 16   0000000000000000000000000000000100011011000000000001101100011011

    TKmerID kmerID = ONE_LSHIFT_63 + dna_encoder("ACGTAAAAACGTACGT");
    TKmerID kmerID_ref = kmerID - ONE_LSHIFT_63;
    std::cout << "ACGTAAAAACGTACGT to bits with head bit:\t" << bits2str(kmerID) << std::endl;
    // head should be 0 after return
    TKmerID kmerID_res = filter_repeats_runs2(kmerID);
    std::cout << "ACGTAAAAACGTACGT to bits after head rem:\t" << bits2str(kmerID) << std::endl;

    if (kmerID_ref != kmerID_res)
    {
        std::cout << "ERROR: expect 0 header, got " << bits2str((MASK_SELECTOR & kmerID_res) >> 54) << std::endl;
        exit(0);
    }
    // |TTTTTCGTAAAAGACGTACGT| = 21
    kmerID = (ONE_LSHIFT_63 >> 5) + dna_encoder("TTTTTCGTAAAAGACGTACGT");
    kmerID_ref = kmerID - (ONE_LSHIFT_63 >> 5);
    kmerID_res = filter_repeats_runs2(kmerID);
    if (kmerID_ref != kmerID_res)
    {
        std::cout << "ERROR: expect zero header, got: " << bits2str((MASK_SELECTOR & kmerID_res) >> 54) << std::endl;
        exit(0);
    }
    std::cout << "INFO: Success, 5Xs filtered\n";

    // case: multiple lengths, only largest filtered out |ACGTACGTACGTAAAA| = 16, |ACGTACGTACGTAAAAATTTT| = 21
    kmerID = ONE_LSHIFT_63 + (ONE_LSHIFT_63 >> 5) + dna_encoder("ACGTACGTACGTAAAAATTTT");
    kmerID_ref = ONE_LSHIFT_63 + dna_encoder("ACGTACGTACGTAAAAATTTT");
    kmerID_res = filter_repeats_runs2(kmerID);
    if (kmerID_res != kmerID_ref)
    {
        std::cout << "ERROR: kmerID_f is 1000000000 ? " << bits2str((MASK_SELECTOR & kmerID_res) >> 54);
        std::cout << " or tail different\n";
    }
    else
        std::cout << "INFO: Success, 5As filtered\n";

    // case: multiple lengths, only largest filtered out |ACGTACGTACGTATATATAT| = 20, |ACGTACGTACGTATATATATA| = 21
    kmerID = (ONE_LSHIFT_63 >> 4) + (ONE_LSHIFT_63 >> 5) + dna_encoder("ACGTACGTACGTATATATATA");
    kmerID_ref = kmerID - (ONE_LSHIFT_63 >> 5);
    //std::cout << "kmerID in: " << std::bitset<64>(kmerID) << " as string " << code2str(kmerID) << std::endl;
    kmerID_res = filter_repeats_runs2(kmerID);

    if (kmerID_ref != kmerID_res)
        std::cout << "ERROR: kmerID_f head should be 0b0000100000, but is: " << bits2str((MASK_SELECTOR & kmerID_res) >> 54) << std::endl;
    else
        std::cout << "INFO: Success, 5TAs filtered\n";

    kmerID = ONE_LSHIFT_63 + (ONE_LSHIFT_63 >> 7) + dna_encoder("CGCGCGCGCGACGTACGTACGTACGT"); // 26
    kmerID_ref = (kmerID - ONE_LSHIFT_63 - (ONE_LSHIFT_63 >> 7)) >> 6; // result trimmed to initially encoded length
    kmerID_res = filter_repeats_runs2(kmerID);
    if (kmerID_ref != kmerID_res)
        std::cout << "ERROR: kmerID_f head should be 0b0000100000, but is: " << bits2str((MASK_SELECTOR & kmerID_res) >> 54) << std::endl;
    else
        std::cout << "INFO: Success, 5TAs filtered\n";

}

int main()
{
    //test_dTM();
    test_filter_repeats_runs();
    return 0;
}
