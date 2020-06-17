#include <bitset>
#include <experimental/filesystem>
#include <iostream>
#include <limits>
#include <vector>

#include "../../src/dna.hpp"
#include "../../src/types/CombinePattern.hpp"
#include "../../src/types/IOConfig.hpp"
#include "../../src/utilities.hpp"

#include "gtest/gtest.h"

namespace fs = std::experimental::filesystem;

using namespace priset;

class io_config_test_f : public ::testing::Test {

private:
    std::string fp = __FILE__;
    std::string rp = fp.substr(0, fp.find_last_of('/'));
    fs::path lib_dir = fs::path(rp) / "../library/3041";
    fs::path work_dir = fs::path(rp) / "../work/3041";

public:
    IOConfig io_cfg;

protected:
    void SetUp() override {
        io_cfg.assign(lib_dir, work_dir, 0, 0);
    }
};


TEST(io_config_test, constructor)
{
    IOConfig io_cfg1;
    IOConfig io_cfg2();
    IOConfig io_cfg3{};
    IOConfig io_cfg4{io_cfg3};
    IOConfig io_cfg5 = io_cfg1;
}

TEST(io_config_test, assign)
{
    std::string fp = __FILE__;
    std::string rp = fp.substr(0, fp.find_last_of('/'));
    fs::path lib_dir = fs::path(rp) / "../library/3041";
    fs::path work_dir = fs::path(rp) / "../work/3041";
    io_cfg.assign(lib_dir, work_dir, 0, 0);
}

TEST_F(io_config_test_f, build_species_set)
{
    EXPECT_EQ(, io_cfg.get_species_count());
}

TEST_F(io_config_test_f, build_taxid2accs_map)
{
    EXPECT_EQ(, io_cfg.get_species_count());
}

TEST_F(io_config_test_f, build_seqNo2acc_map)
{
    EXPECT_EQ(, io_cfg.get_species_count());
}

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
