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
    uint8_t kappa_min;
    uint8_t kappa_max;
    // Tm range default is [52:58]
    std::array<uint8_t, 10> Tm_refs = {44, 46, 50, 54, 56, 58, 60, 64, 68, 70};
    uint64_t prefix_ref = 0b0001110000ULL << 54;
    TKmerID kmerID = (0b1111111111ULL << 54) | dna_encoder("AACGTAACGTAACGTAACGTAACGT");


protected:
    void SetUp() override {
        primer_cfg = new PrimerConfig();
        Tm_min = primer_cfg->get_Tm_min();
        Tm_max = primer_cfg->get_Tm_max();
        kappa_min = primer_cfg->get_kappa_min();
        kappa_max = primer_cfg->get_kappa_max();
    }
};

TEST_F(Tm_test_f, Tm)
{
    uint64_t mask = 1ULL << 63;
    for (auto Tm_ref : Tm_refs)
    {
        EXPECT_EQ(Tm_ref, Tm(kmerID, mask));
        mask >>= 1;
    }
}

TEST_F(Tm_test_f, dTm)
{
    uint64_t mask1 = 1ULL << 63;
    for (uint8_t i = 0; i < Tm_refs.size(); ++i)
    {
        uint64_t mask2 = mask1;
        for (uint8_t j = i; j < Tm_refs.size(); ++j)
        {
            uint8_t dTm_ref = std::abs(Tm_refs[i] - Tm_refs[j]);
            EXPECT_EQ(dTm_ref, dTm(kmerID, mask1, kmerID, mask2));
            mask2 >>= 1;
        }
        mask1 >>= 1;
    }
}

TEST_F(Tm_test_f, Tm_filter)
{
    Tm_filter(kmerID, Tm_min, Tm_max, kappa_min, kappa_max);
    uint64_t prefix = kmerID & PREFIX_SELECTOR;
    EXPECT_EQ(prefix_ref, prefix);
}


int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
