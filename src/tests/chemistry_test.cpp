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
    std::cout << "dTm(" << dna_decoder(kmerID1, mask1) << ", " << dna_decoder(kmerID2, mask2) << ") = " << dTm(kmerID1, mask1, kmerID2, mask2) << std::endl;
    float d = dTm(kmerID1, mask1, kmerID2, mask2);
    if (d != 8)
        std::cout << "ERROR: expect 8, got " << d << std::endl;

    mask1 >>= 2;
    std::cout << "dTm(" << dna_decoder(kmerID1, mask1) << ", " << dna_decoder(kmerID2, mask2) << ") = " << dTm(kmerID1, mask1, kmerID2, mask2) << std::endl;
    mask1 >>= 2;
    std::cout << "dTm(" << dna_decoder(kmerID1, mask1) << ", " << dna_decoder(kmerID2, mask2) << ") = " << dTm(kmerID1, mask1, kmerID2, mask2) << std::endl;
}

void test_filter_repeats_runs()
{
/*    // |ACGTAAAAACGTACGT| = 16
    TKmerID kmerID = ONE_LSHIFT_63 + dna_encoder("ACGTAAAAACGTACGT");
    // head should be 0 after return
    TKmerID kmerID_f = filter_repeats_runs2(kmerID);
    if (!(mask_selector & kmerID_f))
        std::cout << "ERROR: kmerID_f head unequal zero: " << (mask_selector & kmerID_f) << std::endl;
    else
        std::cout << "INFO: Success, 5As filtered\n";
    // |TTTTTCGTAAAAGACGTACGT| = 21
    kmerID = (ONE_LSHIFT_63 >> 5) + dna_encoder("TTTTTCGTAAAAGACGTACGT");
    if (!(mask_selector & (filter_repeats_runs2(kmerID))))
        std::cout << "ERROR: kmerID_f head unequal zero: " << (mask_selector & kmerID_f) << std::endl;
    else
        std::cout << "INFO: Success, 5Ts filtered\n";
*/
    // case: multiple lengths, only largest filtered out |ACGTACGTACGTAAAA| = 16, |ACGTACGTACGTAAAAATTTT| = 21
    TKmerID kmerID = ONE_LSHIFT_63 + (ONE_LSHIFT_63 >> 5) + dna_encoder("ACGTACGTACGTAAAAATTTT");

    //uint64_t kmerID_f = filter_repeats_runs2(kmerID);


    /*if (kmerID_f != kmerID_s)
    {
        std::cout << "ERROR: kmerID_f is 1000000000 ? " << bits2str((mask_selector & kmerID_f) >> 54);
        std::cout << " or tail different\n";
    }
    else
        std::cout << "INFO: Success, 5As filtered\n";
*/
    // case: multiple lengths, only largest filtered out |ACGTACGTACGTATATATAT| = 20, |ACGTACGTACGTATATATATA| = 21
                                                                     //0  ACGTACGTACGTATATATATA
    kmerID = ((0b0000100000ULL + 0b0000010000ULL) << 54) + dna_encoder("ACGTACGTACGTATATATATA");
    //       0000110000000000000001001100110011001100111001001110010011100100
    //09876543210987654321098765432109876543210987654321098765432109876543210
    std::cout << "kmerID in: " << std::bitset<64>(kmerID) << " as string " << code2str(kmerID) << std::endl;
    uint64_t kmerID_f2 = filter_repeats_runs2(kmerID);
    if (!((MASK_SELECTOR & kmerID_f2) & (0b0000100000ULL << 54)))
        std::cout << "ERROR: kmerID_f head should be 0b0000100000, but is: " << bits2str((MASK_SELECTOR & kmerID_f2) >> 54) << std::endl;
    else
        std::cout << "INFO: Success, 5TAs filtered\n";

}

int main()
{
    //test_dTM();
    test_filter_repeats_runs();
    return 0;
}
