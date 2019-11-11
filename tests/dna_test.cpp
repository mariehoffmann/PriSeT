#include "../src/dna.hpp"
#include "../src/primer_cfg_type.hpp"
#include "../src/utilities.hpp"

#include "gtest/gtest.h"

namespace fs = std::experimental::filesystem;
using namespace priset;

void reverse_test()
{
    uint64_t code =     dna_encoder("AGCTAGATGTACTTGT") + ONE_LSHIFT_63;
    uint64_t code_r = dna_encoder("TGTTCATGTAGATCGA") + ONE_LSHIFT_63;
    uint64_t res = reverse(code);
    EXPECT_EQ(complement(code), code_r);
}

void reverse_complement_test()
{
    uint64_t code =     dna_encoder("AGCTAGATGTACTTGT") + ONE_LSHIFT_63;
    uint64_t code_rc = dna_encoder("ACAAGTACATCTAGCT") + ONE_LSHIFT_63;
    auto res = reverse_complement(code);
    EXPECT_EQ(complement(code), code_rc);
}

TEST(dna_test, complement) {
    uint64_t code =   dna_encoder("AGCTAGATGTACTTGT") + ONE_LSHIFT_63;
    uint64_t code_c = dna_encoder("TCGATCTACATGAACA") + ONE_LSHIFT_63;
    uint64_t res = complement(code);
    EXPECT_EQ(complement(code), code_c);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
