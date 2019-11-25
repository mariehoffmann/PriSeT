#include "../src/dna.hpp"
#include "../src/primer_cfg_type.hpp"
#include "../src/utilities.hpp"

#include "gtest/gtest.h"

// -lgtest -lpthread -I /path/to/seqan/include -std=c++14 -O3 -DNDEBUG -W -Wall -pedantic
// g++ ../PriSeT/tests/dna_test.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -lstdc++fs -DNDEBUG -O3 -I/opt/local/include -I/Users/troja/include -L/Users/troja/lib -lsdsl -ldivsufsort -I .. -lgtest -lpthread -o dna_test

namespace fs = std::experimental::filesystem;
using namespace priset;

TEST(dna_test, reverse)
{
    uint64_t code =   dna_encoder("AGCTAGATGTACTTGT") + ONE_LSHIFT_63;
    uint64_t code_r = dna_encoder("TGTTCATGTAGATCGA") + ONE_LSHIFT_63;
    uint64_t res = reverse(code);
    EXPECT_EQ(complement(code), code_r);
}

TEST(dna_test, reverse_complement)
{
    uint64_t code =    dna_encoder("AGCTAGATGTACTTGT") + ONE_LSHIFT_63;
    uint64_t code_rc = dna_encoder("ACAAGTACATCTAGCT") + ONE_LSHIFT_63;
    auto res = reverse_complement(code);
    EXPECT_EQ(complement(code), code_rc);
}

TEST(dna_test, complement)
{
    uint64_t code =   dna_encoder("AGCTAGATGTACTTGT") + ONE_LSHIFT_63;
    uint64_t code_c = dna_encoder("TCGATCTACATGAACA") + ONE_LSHIFT_63;
    uint64_t res = complement(code);
    EXPECT_EQ(complement(code), code_c);
}

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
