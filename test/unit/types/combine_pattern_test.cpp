#include <bitset>
#include <vector>

#include "types/CombinePattern.hpp"
#include "utilities.hpp"

#include "gtest/gtest.h"

using namespace priset;

TEST(combine_pattern_test, constructor)
{
    CombinePattern cp1;
    CombinePattern cp2();
    CombinePattern cp3{};
}

// tests is_set() and none()
TEST(combine_pattern_test, is_set)
{
    CombinePattern cp;
    EXPECT_FALSE(cp.is_set());
    EXPECT_TRUE(cp.none());

    cp.set(1ULL << 63, 1ULL << 63);
    EXPECT_TRUE(cp.is_set());
    EXPECT_FALSE(cp.none());
}

TEST(combine_pattern_test, set)
{
    using pair_type = std::pair<CombinePattern::TOffset, CombinePattern::TOffset>;
    CombinePattern cp;
    cp.set(1ULL << 63, 1ULL << 63);
    cp.set(1ULL << 60, 1ULL << 59);
    cp.set(1ULL << 54, 1ULL << 54);
    std::vector<pair_type> cmb;
    cp.get_combinations(cmb);
    ASSERT_EQ(3ULL, cmb.size());
    pair_type cmb_0{0, 0};
    pair_type cmb_1{3, 4};
    pair_type cmb_2{9, 9};
    EXPECT_EQ(cmb_0.first, cmb.at(0).first);
    EXPECT_EQ(cmb_0.second, cmb.at(0).second);
    EXPECT_EQ(cmb_1.first, cmb.at(1).first);
    EXPECT_EQ(cmb_1.second, cmb.at(1).second);
    EXPECT_EQ(cmb_2.first, cmb.at(2).first);
    EXPECT_EQ(cmb_2.second, cmb.at(2).second);
}

TEST(combine_pattern_test, reset)
{
    using pair_type = std::pair<CombinePattern::TOffset, CombinePattern::TOffset>;
    CombinePattern cp;
    cp.set(1ULL << 63, 1ULL << 63);
    EXPECT_TRUE(cp.is_set());
    cp.set(1ULL << 60, 1ULL << 59);
    EXPECT_TRUE(cp.is_set());

    std::vector<pair_type> cmb;
    pair_type cmb_0{3, 4};

    // remove first combination
    cp.reset(1ULL << 63, 1ULL << 63);
    cp.get_combinations(cmb);
    EXPECT_TRUE(cp.is_set());
    ASSERT_EQ(1ULL, cmb.size());
    EXPECT_EQ(cmb_0, cmb.at(0));

    // remove second combination
    cp.reset(1ULL << 60, 1ULL << 59);
    cp.get_combinations(cmb);
    EXPECT_FALSE(cp.is_set());
    EXPECT_EQ(0ULL, cmb.size());
}

TEST(combine_pattern_test, size)
{
    CombinePattern cp;
    EXPECT_EQ(0ULL, cp.size());
    size_t ctr{0};
    for (uint8_t i = 0; i < 10; ++i)
    {
        for (uint8_t j = 0; j < 10; ++j)
        {
            cp.set(1ULL << (63 - i), 1ULL << (63 - j));
            EXPECT_EQ(++ctr, cp.size());
        }
    }
}

// tests get_combinations and operator[]
TEST(combine_pattern_test, get_combinations)
{
    using pair_type = std::pair<CombinePattern::TOffset, CombinePattern::TOffset>;
    CombinePattern cp;
    size_t ctr{0};
    std::vector<pair_type> cmb;
    for (uint8_t i = 0; i < 10; ++i)
    {
        for (uint8_t j = 0; j < 10; ++j)
        {
            cp.set(1ULL << (63 - i), 1ULL << (63 - j));
            cp.get_combinations(cmb);
            EXPECT_EQ(++ctr, cmb.size());
            pair_type cmb_last{i, j};
            EXPECT_EQ(cmb_last, cmb.back());
            // expect the first i*10 + j bits set to one, rest zero
            for (uint8_t k = 0; k <= i*10 + j; ++k)
                EXPECT_TRUE(cp[k]);
            for (uint8_t k =  i*10 + j + 1; k < 10*10; ++k)
                EXPECT_FALSE(cp[k]);
        }
    }
}

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
