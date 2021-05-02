// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================
//          Author: Marie Hoffmann <ozymandiaz147 AT gmail.com>
//          Manual: https://github.com/mariehoffmann/PriSeT

#pragma once

#include <numeric>
#include <vector>

#include "IOConfig.hpp"
#include "PrimerPair.hpp"
#include "PrimerPairUnpacked.hpp"
#include "PrimerConfig.hpp"

namespace priset
{

template<typename TSeqNoMap>
struct Group
{
    // Default constructor.
    Group() {}

    Group(IOConfig * _io_cfg, TSeqNoMap const * _seqNoMap) :
        io_cfg(_io_cfg), seqNoMap(_seqNoMap){}

    Group(IOConfig * _io_cfg, TSeqNoMap const * _seqNoMap,
        PrimerPairUnpacked<TSeqNoMap> & pair_unpacked) : io_cfg(_io_cfg), seqNoMap(_seqNoMap)
    {
        group_container.push_back(pair_unpacked);
    }

    Group(IOConfig * _io_cfg, TSeqNoMap const * _seqNoMap,
        std::vector<PrimerPairUnpacked<TSeqNoMap>> & _pairs_unpacked) : 
            io_cfg(_io_cfg), 
            seqNoMap(_seqNoMap), 
            group_container(_pairs_unpacked){}

    // Default copy constructor.
    Group(Group const & rhs) = default;

    // Default move constructor.
    Group(Group &&) = default;

    // Copy assignment (non copy-and-swap idiom).
    // TODO: causes error when used in lambda in std::transform
    Group & operator=(Group const & rhs)
    {
        if (this != &rhs)
        {
            io_cfg = rhs.io_cfg;
            seqNoMap = rhs.seqNoMap;
            group_container = rhs.group_container;
            seqNo_cx_vector = rhs.seqNo_cx_vector;
            taxa_count = rhs.taxa_count;
            sequence_count = rhs.sequence_count;
        }
        return *this;
    }

    // Add new pair to group_container.
    void insert(PrimerPairUnpacked<TSeqNoMap> const & pair_unpacked)
    {
        group_container.push_back(pair_unpacked);
        for (TSeqNo seqNo_cx = 0; seqNo_cx < seqNo_cx_vector.size(); ++seqNo_cx)
            seqNo_cx_vector[seqNo_cx] = seqNo_cx_vector.at(seqNo_cx) || pair_unpacked.get_sequence_match(seqNo_cx);
        sequence_count = std::accumulate(seqNo_cx_vector.begin(), seqNo_cx_vector.end(), 0);
    }

    // Return number of taxa amplified by this group of primer pairs.
    size_t get_taxa_count()
    {
        assert(false); // TODO: enable upon provisioning of taxid file
        if (taxa_count)
            return taxa_count;
        std::unordered_set<Taxid> taxa;
        for (uint64_t seqNo_cx = 0; seqNo_cx < seqNo_cx_vector.size(); ++seqNo_cx)
        {
            if (!seqNo_cx_vector.at(seqNo_cx))
                continue;
            TSeqNo seqNo = seqNoMap->at((1ULL << 63) | seqNo_cx);
            Taxid taxid = io_cfg->get_taxid_by_seqNo(seqNo);
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

    std::string get_header()
    {
        return  "#primer forward,primer reverse,frequency,coverage\n";
    }

    std::string to_string() const noexcept
    {
        std::string s = "";
        for (PrimerPairUnpacked<TSeqNoMap> p : group_container)
        {
            s += p.get_forward_primer() + "," + p.get_reverse_primer() + ",";
            s += std::to_string(p.get_frequency()) + "," + std::to_string(p.get_species_count()) + "\n";
        }
        return s;
    }

private:

    // Reference to the I/O configurator.
    IOConfig * io_cfg = nullptr;
    // bidirectional map for reference identifiers as found in the library and
    // continuous idententifiers as used in e.g. in reference lists.
    // seqNoMap : seqNo -> seqNo_cx
    // seqNoMap : 1 << 63 | seqNo_cx -> seqNo
    TSeqNoMap const * seqNoMap = nullptr;

    // PrimerPrimerPairs that are member of this group_container.
    std::vector<PrimerPairUnpacked<TSeqNoMap>> group_container;

    // Bit vector of compressed sequence identifiers.
    std::vector<bool> seqNo_cx_vector;

    // Memoize taxa count.
    size_t taxa_count{0};

    // Memoize reference count
    size_t sequence_count{0};

};

}
