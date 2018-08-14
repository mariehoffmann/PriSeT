// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================

#include <gtest/gtest.h>
#include <priset/types/dna.hpp>

using namespace priset;

TEST(int_types_test, min_viable_uint_v)
{
    auto bool_1_v = detail::min_viable_uint_v<0ull>;
    auto bool_2_v = detail::min_viable_uint_v<1ull>;
    auto uint8_1_v = detail::min_viable_uint_v<2ull>;
    auto uint8_2_v = detail::min_viable_uint_v<0xFFull>;
    auto uint64_2_v = detail::min_viable_uint_v<0xFFFFFFFFFFFFFFFFull>;

    EXPECT_EQ(static_cast<uint64_t>(bool_1_v), 0ull);
    EXPECT_EQ(static_cast<uint64_t>(bool_2_v), 1ull);
    EXPECT_EQ(static_cast<uint64_t>(uint8_1_v), 2ull);
    EXPECT_TRUE((std::is_same_v<decltype(uint16_2_v), uint16_t>));
    EXPECT_TRUE((std::is_same_v<decltype(uint32_1_v), uint32_t>));
    EXPECT_TRUE((std::is_same_v<decltype(uint32_2_v), uint32_t>));
    EXPECT_TRUE((std::is_same_v<decltype(uint64_1_v), uint64_t>));
    EXPECT_TRUE((std::is_same_v<decltype(uint64_2_v), uint64_t>));
}
