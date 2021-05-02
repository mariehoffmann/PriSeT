// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================
//          Author: Marie Hoffmann <ozymandiaz147 AT gmail.com>
//          Manual: https://github.com/mariehoffmann/PriSeT

#pragma once

#include <numeric>
#include <string>
#include <vector>

#include "../dna.hpp"
#include "simple_types.hpp"

namespace priset
{

std::string dna_decoder(uint64_t code);

template<typename TSeqNoMap>
struct PrimerPairUnpacked
{
    PrimerPairUnpacked() = default;
    // Constructor with arguments. Length corrected codes must not be zero.
    PrimerPairUnpacked(IOConfig * _io_cfg, TSeqNoMap const * _seqNo_map,
        std::vector<bool> & _seqNo_cx_vector, uint64_t _code_fwd, uint64_t _code_rev) : 
            io_cfg(_io_cfg), 
            seqNo_map(_seqNo_map), 
            seqNo_cx_vector(_seqNo_cx_vector), 
            code_fwd(_code_fwd), 
            code_rev(_code_rev)
    {
        assert(code_fwd && code_rev);
        for (TSeqNo seqNo_cx = 0; seqNo_cx < seqNo_cx_vector.size(); ++seqNo_cx)
        {
            if (!seqNo_cx_vector.at(seqNo_cx))
                continue;
            // recover original sequence identifier
            // TODO: optionally taxid file can be provided and primer pair information be augmented
            // TSeqNo seqNo = seqNo_map->at((1ULL << 63) | seqNo_cx);
            // species_set.insert(io_cfg->get_taxid_by_seqNo(seqNo));
        }
        it_species = species_set.cbegin();
    }

    // Set flag in reference vector if this pair has a match in a reference
    // with ID = seqNo_cx.
    void set_sequence_match(uint64_t const seqNo_cx)
    {
        // TODO: enable upon provisioning of taxid file
        TSeqNo seqNo = seqNo_map->at((1ULL << 63) | seqNo_cx);
        // TODO: this only upon provisioning of taxid file
        if (false)
        {
            Taxid taxid = io_cfg->get_taxid_by_seqNo(seqNo);
            assert(io_cfg->is_species(taxid));
            species_set.insert(taxid);
        }
        if (seqNo_cx_vector.size() <= seqNo_cx)
            seqNo_cx_vector.resize(seqNo_cx + 1, 0);
        if (!seqNo_cx_vector[seqNo_cx])
            ++seqNo_count;
        seqNo_cx_vector[seqNo_cx] = 1;
        // reset taxid iterator
        it_species = species_set.cbegin();
    }

    // Get flag for match of reference sequence with ID = seqNo_cx.
    bool get_sequence_match(size_t seqNo_cx) const noexcept
    {
        return (seqNo_cx_vector.size() <= seqNo_cx) ? 0 : seqNo_cx_vector[seqNo_cx];
    }

    // Return number of distinct occurences in library.
    size_t get_frequency() noexcept
    {
        if (seqNo_count)
            return seqNo_count;
        seqNo_count = std::accumulate(seqNo_cx_vector.begin(), seqNo_cx_vector.end(), 0);
        return seqNo_count;
    }

    size_t get_species_count() const noexcept
    {
        return species_set.size();
    }

    // Return number of distinct species.
    Taxid get_next_species() noexcept
    {
        if (it_species == species_set.cend())
        {
            it_species = species_set.cbegin();
            return Taxid{0};
        }
        Taxid taxid = *it_species;
        it_species++;
        return taxid;
    }

    // Return vector of compressed sequence identifiers.
    const std::vector<bool> get_seqNo_cx_vector() const noexcept
    {
        return seqNo_cx_vector;
    }

    // Return forward primer as dna string.
    std::string get_forward_primer() const noexcept
    {
        assert(code_fwd);
        return dna_decoder(code_fwd);
    }

    // Return reverse primer as dna string.
    std::string get_reverse_primer() const noexcept
    {
        assert(code_rev);
        return dna_decoder(code_rev);
    }

private:
    // Pointer to I/O configurator.
    IOConfig * io_cfg = nullptr;

    // Pointer to seqNo map resolving sequence IDs and compressed sequence IDs.
    TSeqNoMap const * seqNo_map = nullptr;

    // Bit vector indicating pair occurrences where index corresponds to seqNo_cx.
    std::vector<bool> seqNo_cx_vector;

    // forward code without length information and truncated to true length
    uint64_t code_fwd;

    // reverse code without length information and truncated to true length
    uint64_t code_rev;

    // Memoized sequence counter, i.e. the number of references this pair occurs.
    size_t seqNo_count{0};

    // Set of taxonomic identifiers representing matching references.
    std::unordered_set<Taxid> species_set;

    // Iterator over species_set.
    std::unordered_set<Taxid>::const_iterator it_species;
};

}
