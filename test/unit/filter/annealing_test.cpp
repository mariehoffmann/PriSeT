#include <iostream>

#include "filter/annealing.hpp"
#include "types/all.hpp"
#include "utilities.hpp"

#include "gtest/gtest.h"

using namespace priset;

TEST(self_annealing, connected)
{
    // no connected self-annealing
    uint64_t prefix_ref = 1ULL << 62;
    TKmerID kmerID = prefix_ref | dna_encoder("CGAAACTCAGGGGAACG") ;
    filter_annealing_connected(kmerID);
    uint64_t prefix = kmerID & PREFIX_SELECTOR;
    EXPECT_EQ(prefix_ref, prefix);

    // self-annealing -/- in first 16 positions
    //      CGAAAGTCAGGGGATCG
    //          ||||
    //    CGAAAGTCAGGGGATCG
    kmerID = prefix_ref | dna_encoder("CGAAAGTCAGGGGATCG") ;
    filter_annealing_connected(kmerID);
    prefix = kmerID & PREFIX_SELECTOR;
    EXPECT_EQ(0ULL, prefix);

    // self-annealing -/- in first 16 positions
    //      CGAAAGACAGGGGCTGACT
    //                   ||||
    //    CGAAAGACAGGGGCTGACT
    kmerID = (1ULL << 62) | dna_encoder("CGAAAGTCAGGGGATCG") ;
    filter_annealing_connected(kmerID);
    prefix = kmerID & PREFIX_SELECTOR;
    EXPECT_EQ(0ULL, prefix);

    // self-annealing -/+ below position 16 => discard length bits
    //             TCTAGTCCTCTTCGATCC
    //              ||||
    // CCTAGCTTCTCCTGATCT

    kmerID = (1ULL << 63) | (1ULL << 61) | dna_encoder("TCTAGTCCTCTTCGATCC");
    filter_annealing_connected(kmerID);
    prefix = kmerID & PREFIX_SELECTOR;
    EXPECT_EQ(0ULL, prefix);

    // self-annealing -/+
    //  TGATCGTCTTCGATCCC
    //   |||||    |||||
    // CCCTAGCTTCTGCTAGT
    kmerID = (1ULL << 62) | dna_encoder("TGATCGTCTTCGATCCC");
    filter_annealing_connected(kmerID);
    prefix = kmerID & PREFIX_SELECTOR;
    EXPECT_EQ(0ULL, prefix);

    // reverse self-annealing at positions [14:17] => keep only l = 16
    // 0123456789012345678901234567890
    // TCAGCTCCTCTTCGGATCA
    //               ||||
    //              ACTAGGCTTCTCCTCGACT
    prefix_ref = 1ULL << 63;
    kmerID = (1ULL << 63) | (1ULL << 61) | (1ULL << 60) | dna_encoder("TCAGCACCTCTTCGGATCA");
    filter_annealing_connected(kmerID);
    prefix = kmerID & PREFIX_SELECTOR;
    EXPECT_EQ(prefix_ref, prefix);
}

TEST(self_annealing, disconnected)
{
    // insignifanct self annealing
    //  5-CGAAAGTGAGGGGATCG->
    //    |            | |
    // 5-CGAAAGTGAGGGGATCG->
    uint64_t prefix_ref = 1ULL << 62;
    TKmerID kmerID = prefix_ref | dna_encoder("CGAAAGTGAGGGGATCG");
    filter_annealing_disconnected(kmerID);
    uint64_t prefix = kmerID & PREFIX_SELECTOR;
    EXPECT_EQ(prefix_ref, prefix);

    // annealing positions exceed 50 % when matched -/-
    //    5-AGCTAGATGTACTTGG->
    //      | ||| ||| ||
    // 5-AGCTAGATGTACTTGG->
    kmerID = (1ULL << 63) | dna_encoder("AGCTAGATGTACTTGG");
    filter_annealing_disconnected(kmerID);
    prefix = kmerID & PREFIX_SELECTOR;
    EXPECT_EQ(0ULL, prefix);

    // annealing positions exceed 50 % when matched -/+
    //  CGAAAGACAGGGGCTGACT
    //      || |||   ||| ||
    //      TCAGTCGGGGACAGAAAGC
    kmerID = (1ULL << 63) | (1ULL << 60) | dna_encoder("CGAAAGACAGGGGCTGACT");
    filter_annealing_disconnected(kmerID);
    prefix = kmerID & PREFIX_SELECTOR;
    EXPECT_EQ(0ULL, prefix);

    // annealing positions exceed 50 % when matched -/+
    //   CGAAACTCAGGGGATCGG
    //   |||  | | | |  |||
    //  GGCTAGGGGACTCAAAGC
    kmerID = (1ULL << 61) | dna_encoder("CGAAACTCAGGGGATCGG");
    filter_annealing_disconnected(kmerID);
    prefix = kmerID & PREFIX_SELECTOR;
    EXPECT_EQ(0ULL, prefix);
}

TEST(cross_annealing, connected)
{
    uint64_t prefix, prefix_ref, prefix_ref1, prefix_ref2;
    TKmerID kmerID1, kmerID2;

    // cross-annealing -/- in suffix of 1st code => delete length bit 17 in kmerID1
    // CCGAAAGTCAGGGTGATC  17
    //               ||||
    //             TTCTAGGGCCACGTCT
    prefix_ref = 1ULL << 63;
    kmerID1 = prefix_ref | (1ULL << 62) | dna_encoder("CCGAAAGTCAGGGGATC");
    kmerID2 = prefix_ref | dna_encoder("TTCTAGGGCCACGTCT");
    filter_annealing_connected(kmerID1, kmerID2);
    prefix = kmerID1 & PREFIX_SELECTOR;
    EXPECT_EQ(prefix_ref, prefix);
    prefix = kmerID2 & PREFIX_SELECTOR;
    EXPECT_EQ(prefix_ref, prefix);

    // cross-annealing -/- in both's prefixes => no combination possible
    //         GTAGGATCAGGGGATCG
    //          |||||
    // TTCTAGGGCATCCTCT
    kmerID1 = (1ULL << 63) | (1ULL << 62) |dna_encoder("GTAGGATCAGGGGATCG");
    kmerID2 = (1ULL << 63) | dna_encoder("TTCTAGGGCATCCTCT");
    filter_annealing_connected(kmerID1, kmerID2);
    prefix = kmerID1 & PREFIX_SELECTOR;
    EXPECT_EQ(0ULL, prefix);
    prefix = kmerID2 & PREFIX_SELECTOR;
    EXPECT_EQ(0ULL, prefix);

    // cross-annealing -/+ in both's prefixes
    //       GTAGGATCAGGGGATCG
    //            ||||
    // TCTGCACCGGCTAGTT
    kmerID1 = (1ULL << 63) | (1ULL << 62) | dna_encoder("GTAGGATCAGGGGATCG");
    kmerID2 = (1ULL << 63) | dna_encoder("TTGATCGGCCACGTCT");
    filter_annealing_connected(kmerID1, kmerID2);
    prefix = kmerID1 & PREFIX_SELECTOR;
    EXPECT_EQ(0ULL, prefix);
    prefix = kmerID2 & PREFIX_SELECTOR;
    EXPECT_EQ(0ULL, prefix);

    // cross-annealing -/- between 1st suffix and 2nd prefix => delete length bit 17 in kmerID1
    // cross-annealing -/- between 1st suffix and 2nd suffix => delete length bit 18 in kmerID2
    // 5-GTAGGAACAGAGGATGG->                 5-CTAGGAACAGAGGATGG->
    //                ||||                                 ||||
    //             5-TTACCTGCACCGGCTACT->   5-TTACCTGCACCGGCTACT->
    prefix_ref = 1ULL << 63;
    kmerID1 = prefix_ref | (1ULL << 62) | dna_encoder("CTAGGAACAGAGGATGG");
    kmerID2 = prefix_ref | (1ULL << 61) | dna_encoder("TTACCTGCACCGGCTACT");
    filter_annealing_connected(kmerID1, kmerID2);
    prefix = kmerID1 & PREFIX_SELECTOR;
    EXPECT_EQ(prefix_ref, prefix);
    prefix = kmerID2 & PREFIX_SELECTOR;
    EXPECT_EQ(prefix_ref, prefix);

    // cross-annealing -/- between 1st prefix and 2nd suffix => delete length bit 17 in kmerID2
    //               5-GATCTGCACCGGCTAGTT->
    //                 ||||
    //  5-GTAGGATCAGGGGCTAG->
    //
    // cross-annealing -/+ in both prefixes => delete all length bits
    //          5-GATCTGCACCGGCTAGTT->
    //            ||||
    // <-GATCGGGGACTAGGATG-3
    kmerID1 = (1ULL << 63) | (1ULL << 61) | dna_encoder("GATCTGCACCGGCTAGTT");
    kmerID2 = (1ULL << 63) | (1ULL << 62) | dna_encoder("GTAGGATCAGGGGCTAG");
    filter_annealing_connected(kmerID1, kmerID2);
    prefix = kmerID1 & PREFIX_SELECTOR;
    EXPECT_EQ(0ULL, prefix);
    prefix = kmerID2 & PREFIX_SELECTOR;
    EXPECT_EQ(0ULL, prefix);

    // cross-annealing -/+ between 1st suffix and 2nd prefix => delete length bit 17 in kmerID1
    // 5-GTAGGAACACGGGATCG->
    //                ||||
    // <-GATCTGCACCGGTTAGCC-5
    //
    // cross-annealing +/+ for l'1 = 16 and l'2 = 18 => del length bit 18 of kmerID2
    //   5-CTAGGATCACGGGATCG->
    //                 ||||
    // 5-CCGATTGGCCACGTCTAG->
    prefix_ref = 1ULL << 63;
    kmerID1 = prefix_ref | (1ULL << 62) | dna_encoder("GTAGGAACACGGGATCG");
    kmerID2 = prefix_ref | (1ULL << 61) | dna_encoder("CCGATTGGCCACGTCTAG");
    filter_annealing_connected(kmerID1, kmerID2);
    prefix = kmerID1 & PREFIX_SELECTOR;
    EXPECT_EQ(prefix_ref, prefix);
    prefix = kmerID2 & PREFIX_SELECTOR;
    EXPECT_EQ(prefix_ref, prefix);

    // cross-annealing -/+ between 1st prefix and 2nd suffix => delete length bit 18 in kmerID2
    // <-GATGGGCCATAGAGTAC-5
    //                ||||
    // 5-GTAGGATTCACGGCATGG->
    prefix_ref1 = (1ULL << 63) | (1ULL << 62);
    kmerID1 = prefix_ref1 | dna_encoder("CATGAGATACCGGGTAG");
    prefix_ref2 = 1ULL << 63;
    kmerID2 = prefix_ref2 | (1ULL << 61) | dna_encoder("GTAGGATTCACGGCATGG");
    filter_annealing_connected(kmerID1, kmerID2);
    prefix = kmerID1 & PREFIX_SELECTOR;
    EXPECT_EQ(prefix_ref1, prefix);
    prefix = kmerID2 & PREFIX_SELECTOR;
    EXPECT_EQ(prefix_ref2, prefix);

    // cross-annealing -/- between both suffixes => delete length bit 17 in kmerID1/2
    // note that this decision is arbitrary, more correct would be to generate two
    // solutions - one where only kmerID1 is reset and one where only kmerID2 is reset
    // 5-GATGTGGACCGGCTTGATT->
    //                ||||
    // 5-GTAGGATCACGGAAACT->
    prefix_ref1 = 1ULL << 63;
    kmerID1 = prefix_ref1 | (1ULL << 61) | dna_encoder("GATGTGGACCGGCTTGATT");
    prefix_ref2 = (1ULL << 63) | (1ULL << 62);
    kmerID2 = prefix_ref2 | dna_encoder("GTAGGATCACGGAAACT");
    filter_annealing_connected(kmerID1, kmerID2);
    prefix = kmerID1 & PREFIX_SELECTOR;
    EXPECT_EQ(prefix_ref1, prefix);
    prefix = kmerID2 & PREFIX_SELECTOR;
    EXPECT_EQ(prefix_ref2, prefix);

    // cross-annealing -/+ between both suffixes => delete larger one, i.e. 18 in kmerID1
    // 5-TTGATCGGCCACGTCTAC->
    //                 ||||
    //               <-GATGACCATGCGACCAT-5
    prefix_ref1 = 1ULL << 63;
    kmerID1 = prefix_ref1 | (1ULL << 61) | dna_encoder("TTGATCGGCCACGTCTAC");
    prefix_ref2 = (1ULL << 63) | (1ULL << 62);
    kmerID2 = prefix_ref2 | dna_encoder("TACCAGCGTACCAGTAG");
    filter_annealing_connected(kmerID1, kmerID2);
    prefix = kmerID1 & PREFIX_SELECTOR;
    EXPECT_EQ(prefix_ref1, prefix);
    prefix = kmerID2 & PREFIX_SELECTOR;
    EXPECT_EQ(prefix_ref2, prefix);
}

TEST(cross_annealing, disconnected)
{
    TKmerID kmerID1, kmerID2;
    uint64_t prefix, prefix_ref1, prefix_ref2;

    // insignifanct cross-annealing, not that kmerID2 would not pass self-annealing test!
    //   5-CGAAAGCGAGGGGATCG->
    //       |  |   ||  ||
    // 5-ACATTAACCGGCCGGTAAGGCA->
    prefix_ref1 = 1ULL << 62;
    kmerID1 = prefix_ref1 | dna_encoder("CGAAAGTGAGGGGATCG");
    prefix_ref2 = (1ULL << 62) | (1ULL << 59);
    kmerID2 = prefix_ref2 | dna_encoder("ACATTAACCGGCCGGTAAGGCA");
    filter_annealing_disconnected(kmerID1, kmerID2);
    prefix = kmerID1 & PREFIX_SELECTOR;
    EXPECT_EQ(prefix_ref1, prefix);
    prefix = kmerID2 & PREFIX_SELECTOR;
    EXPECT_EQ(prefix_ref2, prefix);

    // annealing positions exceed 50 % when matched -/-
    //   5-CGAAAGCGAGGGGATCG->
    //       |  ||| ||  ||
    // 5-ACATTAACGCGCCGGTAAGGCA->
    kmerID1 = (1ULL << 62) | dna_encoder("CGAAAGTGAGGGGATCG");
    kmerID2 = (1ULL << 62) | (1ULL << 59) | dna_encoder("ACATTAACGCGCCGGTAAGGCA");
    filter_annealing_disconnected(kmerID1, kmerID2);
    prefix = kmerID1 & PREFIX_SELECTOR;
    EXPECT_EQ(0ULL, prefix);
    prefix = kmerID2 & PREFIX_SELECTOR;
    EXPECT_EQ(0ULL, prefix);


    // annealing positions exceed 50 % when matched -/+
    //  CGAAAGACAGGGGCTGACT
    //       || ||| |  ||
    //      ACTTTCCGCTTCTCTAC
    kmerID1 = (1ULL << 63) | (1ULL << 60) | dna_encoder("CGAAAGACAGGGGCTGACT");
    kmerID2 = (1ULL << 63) | (1ULL << 62) | dna_encoder("CATCTCTTCGCCTTTCA");
    filter_annealing_disconnected(kmerID1, kmerID2);
    prefix = kmerID1 & PREFIX_SELECTOR;
    EXPECT_EQ(0ULL, prefix);
    prefix = kmerID2 & PREFIX_SELECTOR;
    EXPECT_EQ(0ULL, prefix);
}


int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
