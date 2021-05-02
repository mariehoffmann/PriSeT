// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================
//          Author: Marie Hoffmann <ozymandiaz147 AT gmail.com>
//          Manual: https://github.com/mariehoffmann/PriSeT

#pragma once

#include <functional>
#include <sstream>
#include <string>
#include <vector>

#include "../chemistry.hpp"
#include "../dna.hpp"
#include "../filter/CG.hpp"

namespace priset
{

// Result type storing a single primer pair with additional information.
struct Result
{

std::string dna_decoder(uint64_t code, uint64_t const mask);

private:

    // List of distinct taxonomic node identifiers in which this pair occurs.
    std::vector<uint64_t> taxa;

    // List of distinct references in which this pair occurs.
    std::vector<uint64_t> references;

    std::string primer_fwd;
    std::string primer_rev;

    int Tm_fwd{0};
    int Tm_rev{0};

    float CG_fwd{0};
    float CG_rev{0};

public:
    Result(TKmerID kmerID_fwd, uint64_t mask_fwd, TKmerID kmerID_rev, uint64_t mask_rev, std::vector<uint64_t> _taxa, std::vector<uint64_t> _references) :
    taxa(_taxa), references(_references)
    {
        primer_fwd = dna_decoder(kmerID_fwd, mask_fwd);
        primer_rev = dna_decoder(kmerID_rev, mask_rev);
        Tm_fwd = Tm(kmerID_fwd, mask_fwd);
        Tm_rev = Tm(kmerID_rev, mask_rev);
        CG_fwd = CG(kmerID_fwd, mask_fwd);
        CG_rev = CG(kmerID_rev, mask_rev);
    }

    std::pair<std::string, std::string> get_primer_names() const
    {
        std::stringstream ss_fwd, ss_rev;
        ss_fwd << std::hex << std::hash<std::string>{}(primer_fwd);
        ss_rev << std::hex << std::hash<std::string>{}(primer_rev);
        return std::pair<std::string, std::string>{ss_fwd.str(), ss_rev.str()};
    }

    std::pair<std::string, std::string> get_primer_sequences() const
    {
        return std::pair<std::string, std::string>{primer_fwd, primer_rev};
    }

    std::vector<uint64_t> get_taxa() const
    {
        return taxa;
    }

    std::vector<uint64_t> get_references() const
    {
        return references;
    }

    // Return comma-separated result string
    std::string to_string()
    {
        std::stringstream ss;
        ss << primer_fwd << "," << Tm_fwd << "," << CG_fwd << "," << primer_rev << "," << Tm_rev << "," << CG_rev << ",[";
        for (auto taxon : taxa)
        {
            ss << taxon;
            if (taxon != taxa.back())
                ss << ",";
        }
        ss << "]\n";
        return ss.str();
    }
};

}
