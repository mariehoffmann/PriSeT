// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================
//          Author: Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
//          Manual: https://github.com/mariehoffmann/PriSeT

// Primer Settings.

#pragma once

#include <unordered_map>

#include <seqan/basic.h>

//#include "chemistry.hpp"
#include "types.hpp"

namespace priset
{

struct primer_cfg_type
{
    // Member types
    using float_type = float;
    // todo: require member typename TSeq::size_type
    //using TSeq = TSeq; // seqan dna string
    using size_type = uint32_t;

    // interval type for primer length ranges
    using size_interval_type = typename std::pair<size_type, size_type>;

    // The taxonomic identifier type.
    using taxid_type = unsigned int;

    // The index type to identify kmer sequences.
    using kmer_ID_type = unsigned short int;

    // Constructors, destructor and assignment
    // Default constructor.
    constexpr primer_cfg_type() = default;

    // Default copy constructor.
    constexpr primer_cfg_type(primer_cfg_type const &) = default;

    // Default copy construction via assignment.
    primer_cfg_type & operator=(primer_cfg_type const &) = default;

    // Move constructor.
    constexpr primer_cfg_type (primer_cfg_type && rhs) = default;

    // Move assignment.
    primer_cfg_type & operator=(primer_cfg_type && rhs) = default;

    // Use default deconstructor.
    ~primer_cfg_type() = default;
    //!\}

    //
    constexpr void set_root_taxid(taxid_type taxid)
    {
        root_taxid = taxid;
    }

    constexpr taxid_type get_root_taxid() const noexcept
    {
        return root_taxid;
    }

    constexpr float get_occurence_freq() const noexcept
    {
        return occurrence_freq;
    }

    // Set melting temperature range.
    void set_Tm_range(float min_Tm, float max_Tm)
    {
        primer_melt_range = std::make_pair<float, float>(float(min_Tm), float(max_Tm));
    }

    // Set melting temperature range.
    constexpr size_interval_type get_Tm_range() const noexcept
    {
        return primer_melt_range;
    }

    // Get minimal melting temperature.
    constexpr float get_min_Tm() const noexcept
    {
        return primer_melt_range.first;
    }

    // Get maximal melting temperature.
    constexpr float get_max_Tm() const noexcept
    {
        return primer_melt_range.second;
    }

    // Set molar Natrium concentration if salt adjusted method for primer Tm.
    void set_Na(float Na_)
    {
        Na = Na_;
    }

    // Set molar Natrium concentration if salt adjusted method for primer Tm.
    constexpr float get_Na() const noexcept
    {
        return Na;
    }

    // Set primer length range.
    void set_primer_length_range(size_interval_type range_)
    {
        assert(range_.first <= range_.second);
        primer_length_range = range_;
    }


    // Get primer length range
    constexpr size_interval_type get_primer_length_range() const noexcept
    {
        return primer_length_range;
    }

    // Get kmer length, which is the lower bound of the primer length range
    constexpr size_type get_kmer_length() const noexcept
    {
        return primer_length_range.first;
    }

    // Set method for computing primer melting temperature.
    void set_primer_melt_method(enum TMeltMethod method_)
    {
        primer_melt_method = method_;
    }

    // Get method for computing primer melting temperature.
    enum TMeltMethod get_primer_melt_method() const noexcept
    {
        return primer_melt_method;
    }

    // Set distance_range.
    void set_transcript_range(size_interval_type transcript_range_)
    {
        assert(transcript_range_.first >= 50 && transcript_range_.second <= 2000);
        transcript_range = transcript_range_;
    }

    // Get distance_range.
    constexpr size_interval_type get_transcript_range() const noexcept
    {
        return transcript_range;
    }

    // Get CG content range.
    template<typename interval_type>
    void set_CG_content_range(interval_type CG_content_range_) const noexcept
    {
        // assert some boundaries
        CG_content_range = CG_content_range_;
    }

    // Get CG content range.
    template<typename interval_type>
    constexpr interval_type get_CG_content_range() const noexcept
    {
        return CG_content_range;
    }

    // Get lower bound for recommended CG content (relative)
    float get_min_CG_content() const noexcept
    {
        return CG_content_range.first;
    }

    // Get lower bound for recommended CG content (relative)
    float get_max_CG_content() const noexcept
    {
        return CG_content_range.second;
    }

private:

    // Root taxonomic identifier below which references are sampled.
    taxid_type root_taxid{1};

    // Occurrence frequency of primer sequences relative to the number of taxa with at least one existing reference sequence.
    float occurrence_freq{0.1};

    // Default primer melting temperature formular (see chemistry.hpp for details)
    TMeltMethod melt_method{TMeltMethod::WALLACE};

    // Primer length range
    size_interval_type primer_length_range{18, 24};

    // Transcript length range.
    size_interval_type transcript_range{30, 700};

    // Range of primer melting temperatures (best results in range [52-58] degree C).
    std::pair<float, float> primer_melt_range{50.0, 62.0};

    // Maximal permitted temperature difference [Kelvin] of primers.
    // Differences above 5 Kelvin can lead to no amplification.
    float primer_melt_diff{4.0};

    // Method for computing melting temperature of primer.
    TMeltMethod primer_melt_method{TMeltMethod::WALLACE};

    // Molar Natrium concentration for salt-adjusted primer_melt_method.
    float Na{.0};

    // CG content range of primer.
    // It is recommended that the relative content of nts C and G is in the range of 40%-60%.
    std::pair<float, float> CG_content_range{.4, .6};
};

}  // namespace priset
