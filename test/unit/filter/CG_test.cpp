#include "types/PrimerConfig.hpp"
#include "filter/CG.hpp"
#include "utilities.hpp"

#include "gtest/gtest.h"

using namespace priset;

class CG_test_f : public ::testing::Test {

public:
    PrimerConfig * primer_cfg;
    float CG_min;
    float CG_max;
    uint8_t kappa_min;
    uint8_t kappa_max;
    // l = [19,20,24,25] are in range, i.e. CG_content in [.4:.6]
    std::array<float, 10> CG_refs = {6./16., 6./17, 7./18, 8./19, 8./20, 8./21, 8./22, 9./23, 10./24, 10./25};
    TKmerID kmerID = (0b1111111111ULL << 54) | dna_encoder_with_lbit("AACGTAACGTAACGTAACGTAACGT");
    uint64_t prefix_ref = 0b0001100011ULL << 54;


protected:
    void SetUp() override {
        primer_cfg = new PrimerConfig();
        CG_min = primer_cfg->get_CG_min();
        CG_max = primer_cfg->get_CG_max();
        kappa_min = primer_cfg->get_kappa_min();
        kappa_max = primer_cfg->get_kappa_max();
    }
};

TEST_F(CG_test_f, CG_percent)
{
    uint64_t mask = 1ULL << 63;
    for (float CG_ref : CG_refs)
    {
        EXPECT_EQ(CG_ref, CG_percent(kmerID, mask));
        mask >>= 1;
    }
}

TEST_F(CG_test_f, CG_filter)
{
    CG_filter(kmerID, CG_min, CG_max, kappa_min, kappa_max);
    uint64_t prefix = kmerID & PREFIX_SELECTOR;
    EXPECT_EQ(prefix_ref, prefix);
}


int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
