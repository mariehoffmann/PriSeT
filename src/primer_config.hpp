// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================

/*!\file
* \brief Primer Settings. http://www.premierbiosoft.com/tech_notes/PCR_Primer_Design.html
* \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
*/

#pragma once

#include <unordered_map>

#include "chemistry.hpp"
#include "dna.hpp"
#include "interval.hpp"

namespace priset
{

template<typename sequence_type>
struct primer_config
{
    //!\publicsection
    /*!\name Member types
    * \{
    */
    using float_type = float;
    // todo: require member typename sequence_type::size_type
    //using sequence_type = sequence_type; // seqan dna string
    using size_type = typename sequence_type::size_type;
    using taxid_type = unsigned int;
    //!\}

    /* rule of six */
    /*!\name Constructors, destructor and assignment
    * \{
    */
    //!\brief Default constructor.
    constexpr primer_config() = default;

    //!\brief Default copy constructor.
    constexpr primer_config(primer_config const &) = default;

    //!\brief Default copy construction via assignment.
    constexpr primer_config & operator=(primer_config const &) = default;

    //!\brief Move constructor.
    constexpr primer_config (primer_config && rhs) = default;

    //!\brief Move assignment.
    constexpr primer_config & operator=(primer_config && rhs) = default;

    //!\brief Use default deconstructor.
    ~primer_config() = default;
    //!\}

    bool check_primer_constraints(sequence_type p) const noexcept
    {

        return true;
    }

    constexpr bool check_primerpair_constraints(sequence_type p, sequence_type q) const noexcept
    {
        return true;
    }

    //
    constexpr void set_root_taxid(taxid_type taxid)
    {
        root_taxid = taxid;
    }

    constexpr taxid_type get_root_taxid() const noexcept
    {
        return root_taxid;
    }

    //!\brief Set lower bound for percentage of primer-template match.
    constexpr bool set_coverage_rate(float_type coverage_rate_)
    {
        assert(coverage_rate >= .5 && coverage_rate <= 1.0);
        coverage_rate = coverage_rate_;
        return true;
    }

    //!\brief Return lower bound for percentage of primer-template match.
    constexpr float_type get_coverage_rate() const noexcept
    {
        return coverage_rate;
    }

    //!\brief Set mismatch_number.
    constexpr bool set_mismatch_number(size_type mismatch_number_)
    {
        assert(mismatch_number_ >= 0 && mismatch_number_ <= primer_length_range.first);
        mismatch_number = mismatch_number_;
    }

    //!\brief Get mismatch_number.
    constexpr size_type get_mismatch_number() const noexcept
    {
        return mismatch_number;
    }

    //!\brief Set melting temperature range.
    void set_Tm_range(float_type min_Tm, float_type max_Tm)
    {
        primer_melt_range = interval<float_type>{min_Tm, max_Tm};
    }

    //!\brief Set melting temperature range.
    constexpr interval<float_type> get_Tm_range() const noexcept
    {
        return primer_melt_range;
    }

    //!\brief Get minimal melting temperature.
    constexpr float_type get_min_Tm() const noexcept
    {
        return primer_melt_range.min;
    }

    //!\brief Get maximal melting temperature.
    constexpr float_type get_max_Tm() const noexcept
    {
        return primer_melt_range.max;
    }

    //!\brief Set molar Natrium concentration if salt adjusted method for primer Tm.
    bool set_Na(float_type Na_)
    {
        Na = Na_;
        return true;
    }

    //!\brief Set molar Natrium concentration if salt adjusted method for primer Tm.
    constexpr float_type get_Na() const noexcept
    {
        return Na;
    }

    //!\brief Set method for computing primer melting temperature.
    void set_primer_melt_method(chemistry::method method_)
    {
        primer_melt_method = method_;
    }

    //!\brief Get method for computing primer melting temperature.
    constexpr chemistry::method get_primer_melt_method() const noexcept
    {
        return primer_melt_method;
    }

    //!\brief Set distance_range.
    bool set_transcript_range(interval<size_type> transcript_range_)
    {
        assert(transcript_range_.min >= 50 && transcript_range_.max <= 2000);
        transcript_range = transcript_range_;
    }

    //!\brief Get distance_range.
    constexpr interval<size_type> get_transcript_range() const noexcept
    {
        return transcript_range;
    }

    //!\brief Get CG content range.
    void set_CG_content_range(interval<float_type> CG_content_range_) const noexcept
    {
        // assert some boundaries
        CG_content_range = CG_content_range_;
    }

    //!\brief Get CG content range.
    constexpr interval<float_type> get_CG_content_range() const noexcept
    {
        return CG_content_range;
    }

private:

    //!\brief Root taxonomic identifier below which references are sampled.
    taxid_type root_taxid{1};

    //!\brief Lower bound for coverage rate of primer to template alignment.
    float_type coverage_rate{.8};

    //!\brief Upper bound for number of mismatches between primer and template.
    size_type mismatch_number{3};

    //!\brief Primer length range.
    interval<size_type> primer_length_range{interval<size_type>{18, 24}};

    //!\brief Transcript length range.
    interval<size_type> transcript_range{interval<size_type>{30, 700}};

    //!\brief Range of primer melting temperatures (best results in range [52-58] degree C).
    interval<float_type> primer_melt_range{52, 58};

    //!\brief Maximal permitted temperature difference [Kelvin] of primers.
    // Differences above 5 Kelvin can lead to no amplification.
    float_type primer_melt_diff{4.0};

    //!\brief Method for computing melting temperature of primer.
    chemistry::method primer_melt_method{chemistry::method::wallace};

    //!\brief Molar Natrium concentration for salt-adjusted primer_melt_method.
    float_type Na{.0};

    //!\brief CG content range of primer.
    // It is recommended that the relative content of nts C and G is in the range of 40%-60%.
    interval<float_type> CG_content_range{.4, .6};
};

}  // namespace priset
