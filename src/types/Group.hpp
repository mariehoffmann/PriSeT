#pragma once

#include <vector>

#include "Pair.hpp"

namespace priset
{

template<typename TSeqNoMap>
struct Group
{
    Group() = default;

    Group(IOConfig const & _io_cfg, TSeqNoMap const & _seqNoMap) :
    io_cfg(_io_cfg), seqNoMap(_seqNoMap){}

    // Return number of taxa amplified by this group of primer pairs.
    const size_t get_coverage() const constexpr
    {
        if (coverage_count)
            return coverage_count;
        std::unordered_set<Taxid> taxa;
        for (ReferenceID refID : referenceID_set)
        {
            io_cfg.
        }
        return referenceID_set.size();
    }



private:
    IOConfig const & io_cfg;
    // bidirectional map for reference identifiers as found in the library and
    // continuous idententifiers as used in e.g. in reference lists.
    // seqNoMap : seqNo -> seqNo_cx
    // seqNoMap : 1 << 63 | seqNo_cx -> seqNo
    TSeqNoMap const & seqNoMap;

    // Pairs that are member of this group.
    std::vector<PairUnpacked> group;

    // Set of references for which at least on primer pair of this group is effective.
    std::unordered_set<ReferenceID> referenceID_set;

    // Memoize coverage count.
    size_t coverage_count{0};

};

template<typename Pair>
using PairList = std::vector<Pair>;

}
