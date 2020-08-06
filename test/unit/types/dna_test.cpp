#include <bitset>

#include "dna.hpp"
#include "types/PrimerConfig.hpp"
#include "utilities.hpp"

#include "gtest/gtest.h"

using namespace priset;

TEST(dna_test, dna_encoder)
{
    seqan::String<priset::dna> seq = "ACGTACGTACGTACGT";
    uint64_t code_sol = 0b0100011011000110110001101100011011;
    // EXPECT_EQ(code_sol, dna_encoder(seq));
    seq = "AAAACCCCGGGGTTTTGGGGCCCCA";
    code_sol = 0b0100000000010101011010101011111111101010100101010100;
    EXPECT_EQ(code_sol, dna_encoder(seq));
}

TEST(dna_test, dna_encoder_with_lbit)
{
    seqan::String<priset::dna> seq = "ACGTACGTACGTACGT";
    uint64_t code_sol = (1ULL << 63) | 0b0100011011000110110001101100011011;
    EXPECT_EQ(code_sol, dna_encoder_with_lbit(seq));
    seq = "AAAACCCCGGGGTTTTGGGGCCCCA";
    code_sol = (1ULL << 54) | 0b0100000000010101011010101011111111101010100101010100;
    EXPECT_EQ(code_sol, dna_encoder_with_lbit(seq));
}

TEST(dna_test, reverse)
{
    uint64_t code = dna_encoder("AGCTAGATGTACTTGT") + ONE_LSHIFT_63;
    uint64_t code_sol = dna_encoder("TGTTCATGTAGATCGA") + ONE_LSHIFT_63;
    EXPECT_EQ(code_sol, reverse(code));
}

TEST(dna_test, complement)
{
    uint64_t code = dna_encoder("AGCTAGATGTACTTGT") + ONE_LSHIFT_63;
    uint64_t code_sol = dna_encoder("TCGATCTACATGAACA") + ONE_LSHIFT_63;
    EXPECT_EQ(code_sol, complement(code));
}

TEST(dna_test, reverse_complement)
{
    uint64_t code = dna_encoder("AGCTAGATGTACTTGT") + ONE_LSHIFT_63;
    uint64_t code_sol = dna_encoder("ACAAGTACATCTAGCT") + ONE_LSHIFT_63;
    EXPECT_EQ(code_sol, reverse_complement(code));
}

TEST(dna_test, dna_decoder)
{
    // without mask and length = 16
    uint64_t code = (1ULL << 63) | 0b0100011011000110110001101100011011;
    std::string seq_sol = "ACGTACGTACGTACGT";
    EXPECT_EQ(seq_sol, dna_decoder(code));

    // without mask and length = 25
    code = (1ULL << 54) | 0b0100000000010101011010101011111111101010100101010100;
    seq_sol = "AAAACCCCGGGGTTTTGGGGCCCCA";
    EXPECT_EQ(seq_sol, dna_decoder(code));

    // with mask and target length = 16
    code = (1ULL << 63) | 0b0100011011000110110001101100011011;
    seq_sol = "ACGTACGTACGTACGT";
    EXPECT_EQ(seq_sol, dna_decoder(code, 1ULL << 63));

    // with mask and target length = 25
    code = (1ULL << 54) | 0b0100000000010101011010101011111111101010100101010100;
    seq_sol = "AAAACCCCGGGGTTTTGGGGCCCCA";
    EXPECT_EQ(seq_sol, dna_decoder(code, 1ULL << 54));

    // with mask and target length = 19
    code = (1ULL << 54) | 0b0100000000010101011010101011111111101010100101010100;
    seq_sol = "AAAACCCCGGGGTTTTGGG";
    EXPECT_EQ(seq_sol, dna_decoder(code, 1ULL << 60));
}

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
