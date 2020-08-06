#include "types/PrimerConfig.hpp"
#include "filter/Tm.hpp"
#include "utilities.hpp"

#include "gtest/gtest.h"

using namespace priset;

class Tm_test_f : public ::testing::Test {

public:
    PrimerConfig * primer_cfg;
    uint8_t Tm_min;
    uint8_t Tm_max;
    std::array<uint8_t, 10> Tm_refs = {46, 48, 52, 56, 58, 60, 62, 66, 70, 72};
    TKmerID kmerID = (0b1111111111ULL << 54) | dna_encoder("AACGTAACGTAACGTAACGTAACGT");


protected:
    void SetUp() override {
        primer_cfg = new PrimerConfig();
        Tm_min = primer_cfg->get_Tm_min();
        Tm_max = primer_cfg->get_Tm_max();
    }
};

TEST_F(Tm_test_f, Tm)
{
    uint64_t mask = 1ULL << 63;
    for (auto Tm_ref : Tm_refs)
    {
        EXPECT_EQ(Tm_ref, Tm(kmerID, mask));
        mask <<= 1;
    }
}

TEST(Tm_test, dTm)
{
    // kmerID1 = 0001010100|TGCATGCATGCAATGCAATGCAA
    // kmerID2 = 0000010000|AGCATCGATACATCAATCGAT
    TKmerID kmerID1 = 1513342761473819536;
    // TGCATGCATGCAATGCAAT => 11AT | 8CG which is Tm = 54
    uint64_t mask1 = 1ULL << 60;
    TKmerID kmerID2 = 288235407220622179;
    // AGCATCGATACATCAATCGAT => 13AT | 8CG which equals Tm = 58
    uint64_t mask2 = 1ULL << 58;
    std::cout << kmerID2str(kmerID1) << std::endl;
    std::cout << kmerID2str(kmerID2) << std::endl;

    EXPECT_EQ(4, dTm(kmerID1, mask1, kmerID2, mask2));

    // TGCATGCATGCAATGCAATGC => 11AT | 10CG which is Tm = 62
    mask1 >>= 2;
    EXPECT_EQ(4, dTm(kmerID1, mask1, kmerID2, mask2));

    // TGCATGCATGCAATGCAATGCAA => 13AT | 10CG which is Tm = 66
    mask1 >>= 2;
    EXPECT_EQ(8, dTm(kmerID1, mask1, kmerID2, mask2));
}

TEST_F(Tm_test_f, Tm_filter)
{
    // uint64_t prefix_ref = 0b0011100000ULL << 54;
    std::cout << std::bitset<64>(kmerID) << std::endl;
    std::cout << kmerID2str(kmerID) << std::endl;

    Tm_filter(kmerID, primer_cfg->get_Tm_min(), primer_cfg->get_Tm_max());
    uint64_t prefix = kmerID & PREFIX_SELECTOR;
    EXPECT_EQ(prefix_ref, prefix);
}


int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
