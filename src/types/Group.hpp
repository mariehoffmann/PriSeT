#pragma once

#include <numeric>
#include <vector>

#include "IOConfig.hpp"
#include "PrimerPair.hpp"
#include "PrimerConfig.hpp"

namespace priset
{

template<typename TSeqNoMap>
struct Group
{
    // Default constructor.
    Group() {}

    Group(IOConfig const & _io_cfg, TSeqNoMap const & _seqNoMap) :
        io_cfg(_io_cfg), seqNoMap(_seqNoMap){}

    Group(IOConfig const & _io_cfg, TSeqNoMap const & _seqNoMap,
        PrimerPairUnpacked<TSeqNoMap> & pair_unpacked) : io_cfg(_io_cfg), seqNoMap(_seqNoMap)
    {
        group.push_back(pair_unpacked);
    }

    Group(IOConfig const & _io_cfg, TSeqNoMap const & _seqNoMap,
        std::vector<PrimerPairUnpacked<TSeqNoMap>> const & _pairs_unpacked) : io_cfg(_io_cfg), seqNoMap(_seqNoMap)
    {
        group.insert(_pairs_unpacked.cbegin(), _pairs_unpacked.cend(), group.begin());
    }

    // Add new pair to group.
    void insert(PrimerPairUnpacked<TSeqNoMap> const & pair_unpacked)
    {
        group.push_back(pair_unpacked);
        for (TSeqNo seqNo_cx = 0; seqNo_cx < seqNo_cx_vector.size(); ++seqNo_cx)
            seqNo_cx_vector[seqNo_cx] = seqNo_cx_vector.at(seqNo_cx) || pair_unpacked.get_sequence_match(seqNo_cx);
        sequence_count = std::accumulate(seqNo_cx_vector.begin(), seqNo_cx_vector.end(), 0);
    }

    // Return number of taxa amplified by this group of primer pairs.
    size_t get_taxa_count()
    {
        if (taxa_count)
            return taxa_count;
        std::unordered_set<Taxid> taxa;
        for (uint64_t seqNo_cx = 0; seqNo_cx < seqNo_cx_vector.size(); ++seqNo_cx)
        {
            if (!seqNo_cx_vector.at(seqNo_cx))
                continue;
            TSeqNo seqNo = seqNoMap.at((1ULL << 63) | seqNo_cx);
            Accession acc = io_cfg.seqNo2acc_map.at(seqNo);
            Taxid taxid = io_cfg.acc2taxid_map.at(acc);
            taxa.insert(taxid);
        }
        taxa_count = taxa.size();
        return taxa_count;
    }

    // Return number of references amplified by this group of primer pairs.
    size_t get_sequence_count()
    {
        if (sequence_count)
            return sequence_count;
        sequence_count = std::accumulate(seqNo_cx_vector.cbegin(), seqNo_cx_vector.cend(), 0);
        return sequence_count;
    }

    std::string to_string() const noexcept
    {
        std::string s = "#primer forward,primer reverse,frequency,coverage\n";
        for (PrimerPairUnpacked<TSeqNoMap> p : group)
        {
            s += p.get_forward_primer() + "," + p.get_reverse_primer() + ",";
            s += p.get_frequency() + "," + p.get_coverage() + "\n";
        }
        return s;
    }

private:

    // Reference to the I/O configurator.
    IOConfig const & io_cfg;
    // bidirectional map for reference identifiers as found in the library and
    // continuous idententifiers as used in e.g. in reference lists.
    // seqNoMap : seqNo -> seqNo_cx
    // seqNoMap : 1 << 63 | seqNo_cx -> seqNo
    TSeqNoMap const & seqNoMap;

    // PrimerPrimerPairs that are member of this group.
    std::vector<PrimerPairUnpacked<TSeqNoMap>> group;

    // Bit vector of compressed sequence identifiers.
    std::vector<bool> seqNo_cx_vector;

    // Memoize taxa count.
    size_t taxa_count{0};

    // Memoize reference count
    size_t sequence_count{0};

};

}
