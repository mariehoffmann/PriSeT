#include <bitset>
#include <limits>
#include <vector>

#include "../../src/dna.hpp"
#include "../../src/types/CombinePattern.hpp"
#include "../../src/utilities.hpp"

#include "gtest/gtest.h"

using namespace priset;

TEST(primer_pair_test, constructor)
{
    PrimerPair pp1;
    PrimerPair pp2();
    PrimerPair pp3{};
    CombinePattern cp{};
    ASSERT_DEBUG_DEATH((PrimerPair{0, 0, 0, cp}), "");
    PrimerPair pp4{0, 0, 1, cp};
    ASSERT_DEBUG_DEATH((PrimerPair{0, 2, 1, cp}), "");
}

TEST(primer_pair_test, get_seqNo)
{
    PrimerPair pp1;
    ASSERT_DEBUG_DEATH(pp1.get_seqNo(), "");
    CombinePattern cp{};
    PrimerPair pp2{0, 0, 1, cp};
    EXPECT_EQ(0ULL, pp2.get_seqNo());
    PrimerPair pp3{std::numeric_limits<uint64_t>::max(), 2, 3, cp};
    EXPECT_EQ(std::numeric_limits<uint64_t>::max(), pp3.get_seqNo());
}

TEST(primer_pair_test, get_rank_fwd)
{
    PrimerPair pp1;
    ASSERT_DEBUG_DEATH(pp1.get_rank_fwd(), "");
    CombinePattern cp{};
    PrimerPair pp2{0, 0, 1, cp};
    EXPECT_EQ(0ULL, pp2.get_rank_fwd());
    PrimerPair pp3{0, 0, 1, cp};
    EXPECT_EQ(0ULL, pp2.get_rank_fwd());
    PrimerPair pp4{0, 1, std::numeric_limits<uint64_t>::max(), cp};
    EXPECT_EQ(1ULL, pp4.get_rank_fwd());
}

TEST(primer_pair_test, get_rank_rev)
{
    PrimerPair pp1;
    EXPECT_DEBUG_DEATH(pp1.get_rank_rev(), "");
    CombinePattern cp{};
    PrimerPair pp2{0, 0, 1, cp};
    EXPECT_EQ(1ULL, pp2.get_rank_rev());
    PrimerPair pp3{0, 0, std::numeric_limits<uint64_t>::max(), cp};
    EXPECT_EQ(std::numeric_limits<uint64_t>::max(), pp3.get_rank_rev());
}

TEST(primer_pair_test, get_combine_pattern)
{
    CombinePattern cp{};
    cp.set(1ULL << 63, 1ULL << 61);
    PrimerPair pp{0, 0, 1, cp};
    CombinePattern cp_copy = pp.get_combine_pattern();
    EXPECT_EQ(1, cp_copy[2]);
    cp_copy.reset(1ULL << 63, 1ULL << 61);
    EXPECT_TRUE(cp_copy.none());
}

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
