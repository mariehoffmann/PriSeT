

#include "../src/SolverFast.hpp"

#include "gtest/gtest.h"

namespace fs = std::experimental::filesystem;

using namespace priset;

class solver_fast_f : public ::testing::Test {

private:
    using TSeqNoMap = std::unordered_map<TSeqNo, TSeqNo>;
    std::string fp = __FILE__;
    std::string rp = fp.substr(0, fp.find_last_of('/'));
    fs::path lib_dir = fs::path(rp) / "../library/3041";
    fs::path work_dir = fs::path(rp) / "../work/3041";

public:
    IOConfig io_cfg;
    PrimerConfig primer_cfg;
    // TSeqNoMap seqNo_map; // (1 << 63 | seqNo_cx) -> seqNo and seqNo -> seqNo_cx
    // std::vector<bool> seqNo_cx_vector;
    // PrimerPairUnpacked<TSeqNoMap> ppu;
    // std::unordered_set<Taxid> species_set;
    // size_t frequency{0};
    // uint64_t code_fwd = 0b0100000000000000000000000000000000; // A_16
    // uint64_t code_rev = 0b0110101010101010101010101010101010; // G_16

protected:
    void SetUp() override {
        io_cfg.assign(lib_dir, work_dir, 1);
        // set original seqID to compressed 1 -> 0, 3 -> 1, 5 -> 2, 7 -> 3, ...
        // s.t. only every second seqID is registered
        // for (uint32_t seqID_cx = 0; seqID_cx < io_cfg.get_library_size()/2; ++seqID_cx)
        // {
        //     seqNo_map[(seqID_cx << 1) + 1] = (1ULL << 63) | seqID_cx;
        //     seqNo_map[(1ULL << 63) | seqID_cx] = (seqID_cx << 1) + 1;
        // }

        // set every second seqNo_cx to one of first half, i.e. where the pair occurs.
        // seqNo_cx_vector.resize(io_cfg.get_library_size()/4, 0);
        // for (TSeqNo seqNo_cx = 0; seqNo_cx < seqNo_cx_vector.size(); ++++seqNo_cx)
        // {
        //     seqNo_cx_vector.at(seqNo_cx) = 1;
        //     ++frequency;
        //     TSeqNo seqNo = seqNo_map.at((1ULL << 63) | seqNo_cx);
        //     Taxid taxid = io_cfg.get_taxid_by_seqNo(seqNo);
        //     species_set.insert(taxid);
        // }
        // ppu = PrimerPairUnpacked<TSeqNoMap>{&io_cfg, &seqNo_map, seqNo_cx_vector, code_fwd, code_rev};
    }

};

TEST_F(solver_fast_f, constructor)
{
    SolverFast solver = SolverFast(io_cfg, primer_cfg);

}

// TEST_F(solver_fast_f, run)
// {
//     SolverFast solver = SolverFast(io_cfg, primer_cfg);
//     solver.run();
// }
//
// TEST_F(solver_fast_f, group_by_max_coverage)
// {
//     SolverFast solver = SolverFast(io_cfg, primer_cfg);
//     solver.run();
// }
