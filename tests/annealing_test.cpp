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

#include "../src/combine_types.hpp"
#include "../src/filter.hpp"
#include "../src/primer_cfg_type.hpp"
#include "../src/types.hpp"
#include "../src/utilities.hpp"

namespace fs = std::experimental::filesystem;
using namespace priset;

// g++ ../PriSeT/tests/annealing_test.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -lstdc++fs -DNDEBUG -O3 -I/Users/troja/include -L/Users/troja/lib -lsdsl -ldivsufsort -I .. -o annealing_test

void test_filter_self_annealing_connected()
{
    // no connected self-annealing
    TKmerID kmerID = (1ULL << 62) | dna_encoder("CGAAACTCAGGGGAACG") ;
    filter_annealing_connected(kmerID);
    if (!(PREFIX_SELECTOR & kmerID))
        std::cout << "ERROR: expected kmerID to pass annealing test\n";
    else
        std::cout << "OK\n";

    // self-annealing -/- in first 16 positions
    //      CGAAAGTCAGGGGATCG
    //          ||||
    //    CGAAAGTCAGGGGATCG
    kmerID = (1ULL << 62) | dna_encoder("CGAAAGTCAGGGGATCG") ;
    filter_annealing_connected(kmerID);
    if ((PREFIX_SELECTOR & kmerID))
        std::cout << "ERROR: expected kmerID not to pass annealing test\n";
    else
        std::cout << "OK\n";

    // self-annealing -/- in first 16 positions
    //      CGAAAGACAGGGGCTGACT
    //                   ||||
    //    CGAAAGACAGGGGCTGACT
    kmerID = (1ULL << 62) | dna_encoder("CGAAAGTCAGGGGATCG") ;
    filter_annealing_connected(kmerID);
    if ((PREFIX_SELECTOR & kmerID))
        std::cout << "ERROR: expected kmerID not to pass annealing test\nt";
    else
        std::cout << "OK\n";


    // self-annealing -/+ below position 16 => discard length bits
    //             TCTAGTCCTCTTCGATCC
    //              ||||
    // CCTAGCTTCTCCTGATCT
    kmerID = (1ULL << 63) | (1ULL << 61) | dna_encoder("TCTAGTCCTCTTCGATCC");
    filter_annealing_connected(kmerID);
    if (PREFIX_SELECTOR & kmerID)
        std::cout << "ERROR: expected kmerID not to pass annealing test\n";
    else
        std::cout << "OK\n";

    // self-annealing -/+
    //  TGATCGTCTTCGATCCC
    //   |||||    |||||
    // CCCTAGCTTCTGCTAGT
    kmerID = dna_encoder("TGATCGTCTTCGATCCC") + (ONE_LSHIFT_63 >> 2);
    filter_annealing_connected(kmerID);
    if (PREFIX_SELECTOR & kmerID)
        std::cout << "ERROR: expected kmerID not to pass annealing test\n";
    else
        std::cout << "OK\n";

    // reverse self-annealing below position 16 => discard length bits
    // 0123456789012345678901234567890
    // TCAGCTCCTCTTCGGATCA
    //               ||||
    //              ACTAGGCTTCTCCTCGACT
    kmerID = (1ULL << 63) | (1ULL << 61) | (1ULL << 60) | dna_encoder("TCAGCACCTCTTCGGATCA");
    filter_annealing_connected(kmerID);
    std::cout << std::bitset<10>(kmerID >> 53) << std::endl;
    if ((PREFIX_SELECTOR & kmerID) != (1ULL << 63))
        std::cout << "ERROR: expected kmerID to pass annealing test only for length 16\n";
    else
        std::cout << "OK\n";
}

void test_filter_self_annealing_disconnected()
{
    // CGAAAGTCAGGGGATCG
    //           |  |  |
    //           CGATCCCCTGACTTTCG
    TKmerID kmerID = (1ULL << 62) | dna_encoder("CGAAAGTCAGGGGATCG");
    //filter_annealing_connected(kmerID);
    if (!(PREFIX_SELECTOR & kmerID))
        std::cout << "ERROR: expected kmerID to pass annealing test\n";
    else
        std::cout << "OK\n";

    // Reverse self-annealing
    //  CGAAAGACAGGGGCTGACT
    //      || |||   ||| ||
    //      TCAGTCGGGGACAGAAAGC
    kmerID = (1ULL << 63) | (1ULL << 60) | dna_encoder("CGAAAGACAGGGGCTGACT");
    //filter_annealing_connected(kmerID);
    if ((PREFIX_SELECTOR & kmerID) != (1ULL << 63))
        std::cout << "ERROR: expected kmerID not to pass annealing test only for length 19\n";
    else
        std::cout << "OK\n";

    // Reverse self-annealing
    //  CGAAACTCAGGGGATCG
    //  |||  | | | |  |||
    //  GCTAGGGGACTCAAAGC

        /*
        * case: accumulated annealing positions exceed 50 % when matched -/-
           5-AGCTAGATGTACTTG->
             | ||| ||| ||
        5-AGCTAGATGTACTTG->
        */
        // TKmerID kmerID4 = dna_encoder("AGCTAGATGTACTTGT") + ONE_LSHIFT_63;
        // filter_annealing_connected(kmerID4);
        // if (PREFIX_SELECTOR & kmerID4)
        //     std::cout << "ERROR: expected kmerID not to pass annealing test\nt";
        // else
        //     std::cout << "OK\n";

    /* case: accumulated annealing positions exceed 50 % when matched +/-
         5-ACTTAGATGTACGTGG->
           ||  || || || ||
       <-GGTGCATGTAGATTCA-5
    */
    /*    TKmerID kmerID5 = dna_encoder("ACTTAGATGTACGTGG") + (ONE_LSHIFT_63 >> 1);
        filter_annealing_connected(kmerID5);
        if (PREFIX_SELECTOR & kmerID5)
            std::cout << "ERROR: expected kmerID not to pass annealing test\nt";
        else
            std::cout << "OK\n";
    */

}

void test_filter_cross_annealing_connected()
{
    // cross-annealing -/- in suffix of 1st code => delete length bit 17 in kmerID1
    // CGAAAGTCAGGGTGATC  17
    //  TTCTAGGGCCACGTCT
    //              ||||
    //            TTCTAGGGCCACGTCT

    TKmerID kmerID1 = (1ULL << 63) | (1ULL << 62) | dna_encoder("CGAAAGTCAGGGGATC");
    TKmerID kmerID2 = (1ULL << 63) | dna_encoder("TTCTAGGGCCACGTCT");
    filter_annealing_connected(kmerID1, kmerID2);
    if ((PREFIX_SELECTOR & kmerID1) != (1ULL << 63) || (PREFIX_SELECTOR & kmerID2) != (1ULL << 63))
        std::cout << "ERROR: expected kmerIDs to pass annealing test\nt";
    else
        std::cout << "OK\n";

    // case: forward cross annealing in both's prefixes => no combination possible
    //         GTAGGATCAGGGGATCG
    //          |||||
    // TTCTAGGGCATCCTCT
    // kmerID1 = (1ULL << 63) | (1ULL << 62) |dna_encoder("CGAAAGTCAGGGGATCG");
    // kmerID2 = dna_encoder("TTCTAGGGCCACGTCT") + ONE_LSHIFT_63;
    // filter_annealing_connected(kmerID1, kmerID2);
    // if (!(PREFIX_SELECTOR & kmerID1) || !(PREFIX_SELECTOR & kmerID2))
    //     std::cout << "ERROR: expected kmerIDs not to pass annealing test\nt";
    // else
    //     std::cout << "OK\n";

    // // case: reverse cross annealing in both's prefixes
    // kmerID2 = dna_encoder("TTGATCGGCCACGTCT") + ONE_LSHIFT_63;
    // filter_annealing_connected(kmerID1, kmerID2);
    // if (!(PREFIX_SELECTOR & kmerID1) || !(PREFIX_SELECTOR & kmerID2))
    //     std::cout << "ERROR: expected kmerIDs not to pass annealing test\nt";
    // else
    //     std::cout << "OK\n";
    //
    // // case: forward cross annealing between 1st prefix and 2nd suffix
    // kmerID2 = dna_encoder("TTGATCGGCCACGTCTAG") + ONE_LSHIFT_63 + (ONE_LSHIFT_63 >> 2);
    // filter_annealing_connected(kmerID1, kmerID2);
    // if (!(PREFIX_SELECTOR & kmerID1) && !(PREFIX_SELECTOR & kmerID2))
    //     std::cout << "ERROR: expected kmerID2 not to pass annealing test\nt";
    // else
    //     std::cout << "OK\n";

    // case: reverse cross annealing between 1st prefix and 2nd suffix
    // case: forward cross annealing between 1st suffix and 2nd prefix
    // case: reverse cross annealing between 1st suffix and 2nd prefix
    // case: forward annealing both's suffixes
    // case: reverse annealing both's suffixes

}

void test_filter_cross_annealing_disconnected()
{

}

int main()
{
    //test_filter_self_annealing_connected();
     test_filter_cross_annealing_connected();

    return 0;
}
