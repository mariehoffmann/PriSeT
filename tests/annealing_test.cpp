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
    // insignifanct self annealing
    //  5-CGAAAGTGAGGGGATCG->
    //    |            | |
    // 5-CGAAAGTGAGGGGATCG->
    TKmerID kmerID = (1ULL << 62) | dna_encoder("CGAAAGTGAGGGGATCG");
    filter_annealing_disconnected(kmerID);
    if (!(PREFIX_SELECTOR & kmerID))
        std::cout << "ERROR: expected kmerID to pass annealing test\n";
    else
        std::cout << "OK\n";

    // annealing positions exceed 50 % when matched -/-
    //    5-AGCTAGATGTACTTGG->
    //      | ||| ||| ||
    // 5-AGCTAGATGTACTTGG->
    kmerID = dna_encoder("AGCTAGATGTACTTGG") + ONE_LSHIFT_63;
    filter_annealing_disconnected(kmerID);
    if (PREFIX_SELECTOR & kmerID)
        std::cout << "ERROR: expected kmerID not to pass annealing test\nt";
    else
        std::cout << "OK\n";

    // annealing positions exceed 50 % when matched -/+
    //  CGAAAGACAGGGGCTGACT
    //      || |||   ||| ||
    //      TCAGTCGGGGACAGAAAGC
    kmerID = (1ULL << 63) | (1ULL << 60) | dna_encoder("CGAAAGACAGGGGCTGACT");
    filter_annealing_disconnected(kmerID);
    if (PREFIX_SELECTOR & kmerID)
        std::cout << "ERROR: expected kmerID not to pass annealing test\n";
    else
        std::cout << "OK\n";

    // annealing positions exceed 50 % when matched -/+
    //   CGAAACTCAGGGGATCGG
    //   |||  | | | |  |||
    //  GGCTAGGGGACTCAAAGC
    // 00000000001000000000
    kmerID = (1ULL << 61) | dna_encoder("CGAAACTCAGGGGATCGG");
    filter_annealing_disconnected(kmerID);
    if (PREFIX_SELECTOR & kmerID)
        std::cout << "ERROR: expected kmerID not to pass annealing test.\n";
    else
        std::cout << "OK\n";
}

void test_filter_cross_annealing_connected()
{
    // cross-annealing -/- in suffix of 1st code => delete length bit 17 in kmerID1
    // CCGAAAGTCAGGGTGATC  17
    //               ||||
    //             TTCTAGGGCCACGTCT

    TKmerID kmerID1 = (1ULL << 63) | (1ULL << 62) | dna_encoder("CCGAAAGTCAGGGGATC");
    TKmerID kmerID2 = (1ULL << 63) | dna_encoder("TTCTAGGGCCACGTCT");
    filter_annealing_connected(kmerID1, kmerID2);
    if ((PREFIX_SELECTOR & kmerID1) != (1ULL << 63) || (PREFIX_SELECTOR & kmerID2) != (1ULL << 63))
        std::cout << "ERROR: expected kmerIDs to pass annealing test\nt";
    else
        std::cout << "OK\n";

    // // cross-annealing -/- in both's prefixes => no combination possible
    // //         GTAGGATCAGGGGATCG
    // //          |||||
    // // TTCTAGGGCATCCTCT
    kmerID1 = (1ULL << 63) | (1ULL << 62) |dna_encoder("GTAGGATCAGGGGATCG");
    kmerID2 = dna_encoder("TTCTAGGGCATCCTCT") + ONE_LSHIFT_63;
    filter_annealing_connected(kmerID1, kmerID2);
    if ((PREFIX_SELECTOR & kmerID1) || (PREFIX_SELECTOR & kmerID2))
        std::cout << "ERROR: expected kmerIDs not to pass annealing test\nt";
    else
        std::cout << "OK\n";

    // // cross-annealing -/+ in both's prefixes
    //       GTAGGATCAGGGGATCG
    //            ||||
    // TCTGCACCGGCTAGTT
    kmerID1 = (1ULL << 63) | (1ULL << 62) | dna_encoder("GTAGGATCAGGGGATCG");
    kmerID2 = (1ULL << 63) | dna_encoder("TTGATCGGCCACGTCT");
    filter_annealing_connected(kmerID1, kmerID2);
    if ((PREFIX_SELECTOR & kmerID1) || (PREFIX_SELECTOR & kmerID2))
        std::cout << "ERROR: expected kmerIDs not to pass annealing test\nt";
    else
        std::cout << "OK\n";

    // cross-annealing -/- between 1st suffix and 2nd prefix => delete length bit 17 in kmerID1
    // cross-annealing -/- between 1st suffix and 2nd suffix => delete length bit 18 in kmerID2
    // 5-GTAGGAACAGAGGATGG->                 5-CTAGGAACAGAGGATGG->
    //                ||||                                 ||||
    //             5-TTACCTGCACCGGCTACT->   5-TTACCTGCACCGGCTACT->
    kmerID1 = (1ULL << 63) | (1ULL << 62) | dna_encoder("CTAGGAACAGAGGATGG");
    kmerID2 = (1ULL << 63) | (1ULL << 61) | dna_encoder("TTACCTGCACCGGCTACT");
    filter_annealing_connected(kmerID1, kmerID2);
    if (((PREFIX_SELECTOR & kmerID1) != (1ULL << 63)) || (PREFIX_SELECTOR & kmerID2) != (1ULL << 63))
        std::cout << "ERROR: expected kmerID2 not to pass annealing test\nt";
    else
        std::cout << "OK\n";

    // cross-annealing -/- between 1st prefix and 2nd suffix => delete length bit 17 in kmerID2
    //               5-GATCTGCACCGGCTAGTT->
    //                 ||||
    //  5-GTAGGATCAGGGGCTAG->
    kmerID1 = (1ULL << 63) | (1ULL << 61) | dna_encoder("GAACTGCACCGGCTAGTA");
    kmerID2 = (1ULL << 63) | (1ULL << 62) | dna_encoder("GAAGGATCAGCGGCAAG");
    filter_annealing_connected(kmerID1, kmerID2);
    if (((PREFIX_SELECTOR & kmerID1) == ((1ULL << 63) | (1ULL << 61))) && (PREFIX_SELECTOR & kmerID2) == (1ULL << 63))
        std::cout << "ERROR: expected l1 = 17 to be deleted\n";
    else
        std::cout << "OK\n";

    // cross-annealing -/+ between 1st suffix and 2nd prefix => delete length bit 17 in kmerID1
    // 5-GTAGGAACACGGGATCG->
    //                ||||
    // <-GATCTGCACCGGTTAGCC-5
    //
    // cross-annealing +/+ for l'1 = 16 and l'2 = 18 => del length bit 18 of kmerID2
    //           5-CTAGGATCACGGGATCG->
    //                         ||||
    //         5-CCGATTGGCCACGTCTAG->
    kmerID1 = (1ULL << 63) | (1ULL << 62) | dna_encoder("GTAGGAACACGGGATCG");
    kmerID2 = (1ULL << 63) | (1ULL << 61) | dna_encoder("CCGATTGGCCACGTCTAG");
    filter_annealing_connected(kmerID1, kmerID2);
    if (((PREFIX_SELECTOR & kmerID1) != (1ULL << 63)) || (PREFIX_SELECTOR & kmerID2) != (1ULL << 63))
        std::cout << "ERROR: expected kmerID2 not to pass annealing test\nt";
    else
        std::cout << "OK\n";

    // cross-annealing -/+ between 1st prefix and 2nd suffix => delete length bit 18 in kmerID2
    // <-GATGGGCCATAGAGTAC-5
    //                ||||
    // 5-GTAGGATTCACGGCATGG->
    kmerID1 = (1ULL << 63) | (1ULL << 62) | dna_encoder("CATGAGATACCGGGTAG");
    kmerID2 = (1ULL << 63) | (1ULL << 61) | dna_encoder("GTAGGATTCACGGCATGG");
    filter_annealing_connected(kmerID1, kmerID2);
    if (((PREFIX_SELECTOR & kmerID1) != ((1ULL << 63) | (1ULL << 62))) || (PREFIX_SELECTOR & kmerID2) != (1ULL << 63))
        std::cout << "ERROR: expected kmerID2 not to pass annealing test\nt";
    else
        std::cout << "OK\n";

    // cross-annealing -/- between both suffixes => delete length bit 17 in kmerID1/2
    // 5-GATGTGGACCGGCTTGAT->
    //                ||||
    // 5-GTAGGATCACGGAAACT->

//    GATGTGGACCGGCTTGAT
//                  ||||
//          TCAAAGGCACTAGGATG

    kmerID1 = (1ULL << 63) | (1ULL << 61) | dna_encoder("GATGTGGACCGGCTTGAT");
    kmerID2 = (1ULL << 63) | (1ULL << 62) | dna_encoder("GTAGGATCACGGAAACT");
    filter_annealing_connected(kmerID1, kmerID2);
    if (((PREFIX_SELECTOR & kmerID1) != (1ULL << 63)) || (PREFIX_SELECTOR & kmerID2) != ((1ULL << 63) | (1ULL << 62)))
        std::cout << "ERROR: expected kmerID2 not to pass annealing test\nt";
    else
        std::cout << "OK\n";

    // cross-annealing -/+ between both suffixes => delete larger one, i.e. 18 in kmerID1
    // 5-TTGATCGGCCACGTCTAC->
    //                 ||||
    //               <-GATGACCATGCGACCAT-5
    kmerID1 = (1ULL << 63) | (1ULL << 61) | dna_encoder("TTGATCGGCCACGTCTAC");
    kmerID2 = (1ULL << 63) | (1ULL << 62) | dna_encoder("TACCAGCGTACCAGTAG");
    filter_annealing_connected(kmerID1, kmerID2);
    if (((PREFIX_SELECTOR & kmerID1) != ((1ULL << 63))) || (PREFIX_SELECTOR & kmerID2) != ((1ULL << 63) | (1ULL << 62)))
        std::cout << "ERROR: expected kmerID2 not to pass annealing test\nt";
    else
        std::cout << "OK\n";
}

void test_filter_cross_annealing_disconnected()
{

}

int main()
{
    //test_filter_self_annealing_connected();
    //test_filter_cross_annealing_connected();
    test_filter_self_annealing_disconnected();
    return 0;
}
