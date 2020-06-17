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

class primer_pair_unpacked_test_f : public ::testing::Test {

private:
    using TSeqNoMap = std::unordered_map<TSeqNo, TSeqNo>;
    std::string fp = __FILE__;
    std::string rp = fp.substr(0, fp.find_last_of('/'));
    fs::path lib_dir = fs::path(rp) / "../library/3041";
    fs::path work_dir = fs::path(rp) / "../work/3041";

public:
    IOConfig io_cfg;
    TSeqNoMap seqNo_map; // (1 << 63 | seqNo_cx) -> seqNo and seqNo -> seqNo_cx
    std::vector<bool> seqNo_cx_vector;
    PrimerPairUnpacked<TSeqNoMap> ppu;
    std::unordered_set<Taxid> species_set;
    size_t frequency{0};
    uint64_t code_fwd = 0b0100000000000000000000000000000000; // A_16
    uint64_t code_rev = 0b0110101010101010101010101010101010; // G_16

protected:
    void SetUp() override {
        io_cfg.assign(lib_dir, work_dir, 0, 0);
        // set original seqID to compressed 1 -> 0, 3 -> 1, 5 -> 2, 7 -> 3, ...
        // s.t. only every second seqID is registered
        for (uint32_t seqID_cx = 0; seqID_cx < io_cfg.get_library_size()/2; ++seqID_cx)
        {
            seqNo_map[(seqID_cx << 1) + 1] = (1ULL << 63) | seqID_cx;
            seqNo_map[(1ULL << 63) | seqID_cx] = (seqID_cx << 1) + 1;
        }

        // set every second seqNo_cx to one of first half, i.e. where the pair occurs.
        seqNo_cx_vector.resize(io_cfg.get_library_size()/4, 0);
        for (TSeqNo seqNo_cx = 0; seqNo_cx < seqNo_cx_vector.size(); ++++seqNo_cx)
        {
            seqNo_cx_vector.at(seqNo_cx) = 1;
            ++frequency;
            TSeqNo seqNo = seqNo_map.at((1ULL << 63) | seqNo_cx);
            Taxid taxid = io_cfg.get_taxid_by_seqNo(seqNo);
            species_set.insert(taxid);
        }
        ppu = PrimerPairUnpacked<TSeqNoMap>{&io_cfg, &seqNo_map, seqNo_cx_vector, code_fwd, code_rev};
    }

};


TEST_F(primer_pair_unpacked_test_f, constructor)
{
    using TSeqNoMap = std::unordered_map<TSeqNo, TSeqNo>;
    PrimerPairUnpacked<TSeqNoMap> ppu1;
    PrimerPairUnpacked<TSeqNoMap> pp2();
    PrimerPairUnpacked<TSeqNoMap> pp3{};
    // initializer list construction with empty structures
    IOConfig io_cfg;
    TSeqNoMap seqNo_map;
    std::vector<bool> seqNo_cx_vector;
    PrimerPairUnpacked<TSeqNoMap> pp4{&io_cfg, &seqNo_map, seqNo_cx_vector, code_fwd, code_rev};
}

TEST_F(primer_pair_unpacked_test_f, get_seqNo_cx_vector)
{
    std::vector<bool> seqNo_cx_vector_cp = ppu.get_seqNo_cx_vector();
    EXPECT_EQ(seqNo_cx_vector, seqNo_cx_vector_cp);
}

TEST_F(primer_pair_unpacked_test_f, set_sequence_match)
{
    std::vector<bool> seqNo_cx_vector1 = ppu.get_seqNo_cx_vector();
    std::cout << "size of seqNo_cx_vector1 = " << seqNo_cx_vector1.size() << std::endl;

    // set sequence match existing before
    ppu.set_sequence_match(2);
    std::vector<bool> seqNo_cx_vector2 = ppu.get_seqNo_cx_vector();
    EXPECT_EQ(seqNo_cx_vector1, seqNo_cx_vector2);

    // set sequence match not existing before
    EXPECT_EQ(0, seqNo_cx_vector2.at(5));
    ppu.set_sequence_match(5);
    seqNo_cx_vector2 = ppu.get_seqNo_cx_vector();
    EXPECT_EQ(1, seqNo_cx_vector2.at(5));
    seqNo_cx_vector2[5] = 0;
    EXPECT_EQ(seqNo_cx_vector2, seqNo_cx_vector1);

    // force resizing of bit vector
    ppu.set_sequence_match(seqNo_cx_vector1.size() + 1);
    seqNo_cx_vector2 = ppu.get_seqNo_cx_vector();
    EXPECT_EQ(0, seqNo_cx_vector2.at(seqNo_cx_vector1.size()));
    EXPECT_EQ(1, seqNo_cx_vector2.at(seqNo_cx_vector1.size()) + 1);
}

TEST_F(primer_pair_unpacked_test_f, get_sequence_match)
{
    EXPECT_EQ(1, ppu.get_sequence_match(0));
    EXPECT_EQ(0, ppu.get_sequence_match(1));
}

TEST_F(primer_pair_unpacked_test_f, get_frequency)
{
    EXPECT_EQ(frequency, ppu.get_frequency());
    ppu.set_sequence_match(seqNo_cx_vector.size());
    EXPECT_EQ(frequency + 1, ppu.get_frequency());
}

TEST_F(primer_pair_unpacked_test_f, get_coverage)
{
    EXPECT_EQ(species_set.size(), ppu.get_coverage());
}

TEST_F(primer_pair_unpacked_test_f, get_forward_primer)
{
    EXPECT_EQ(std::string(16, 'A'), ppu.get_forward_primer());
}

TEST_F(primer_pair_unpacked_test_f, get_reverse_primer)
{
    EXPECT_EQ(std::string(16, 'G'), ppu.get_reverse_primer());
}

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
