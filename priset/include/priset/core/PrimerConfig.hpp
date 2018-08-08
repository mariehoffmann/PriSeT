// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================

/*!\file
 * \brief Primer Settings.
 * \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 */

 #pragma once

 #include <unordered_map>

 #include <priset/core/chemistry.hpp>
 #include <priset/types/dna.hpp>


// satisfies the primer_config_concept.
namespace priset
{

template<typename value_type>
struct Range
{
    value_type min;
    value_type max;
};

struct PrimerConfig
{
    //!\publicsection
    /*!\name Member types
    * \{
    */
    using float_type = float;
    using sequence_type = std::String; // seqan dna string
    using size_type = unsigned short int;
    using taxid_type = unsigned integer;
    //!\}

    /* rule of six */
    /*!\name Constructors, destructor and assignment
    * \{
    */
    //!\brief Default constructor.
    constexpr PrimerConfig() :
        _root_taxid{1},
        _coverage_rate{.8},
        _mismatch_number{3},
        _transcript_range{Range<size_type>{30, 700}},
        _primer_length_range{Range<size_type>{18, 24}},
        _primer_melt_range{Range<float_type>{50.0, 60.0}},
        _primer_melt_diff{10.0},
        _primer_melt_method{chemistry::wallace},
        _Na{.0},
        {};



    //!\brief Default copy constructor.
    constexpr PrimerConfig(PrimerConfig const &) = default;

    //!\brief Default copy construction via assignment.
    constexpr PrimerConfig & operator=(PrimerConfig const &) = default;

    //!\brief Move constructor.
    constexpr PrimerConfig (PrimerConfig && rhs) = default;

    //!\brief Move assignment.
    constexpr PrimerConfig & operator=(PrimerConfig && rhs) = default;

    //!\brief Use default deconstructor.
    ~PrimerConfig() = default;
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
    constexpr bool set_root_taxid(taxid_type taxid)
    {
        _root_taxid = taxid;
    }

    constexpr taxid_type get_root_taxid() const noexcept
    {
        return _root_taxid;
    }

    //!\brief Set lower bound for percentage of primer-template match.
    constexpr bool set_coverage_rate(float_type coverage_rate)
    {
        assert(coverage_rate >= .5 && coverage_rate <= 1.0);
        _coverage_rate = coverage_rate;
        return true;
    }

    //!\brief Return lower bound for percentage of primer-template match.
    constexpr float_type get_coverage_rate() const noexcept
    {
        return _coverage_rate;
    }

    //!\brief Set mismatch_number.
    constexpr bool set_mismatch_number(size_type mismatch_number)
    {
        assert(mismatch_number >= 0 && mismatch_number <= _primer_length_range.first);
        _mismatch_number = mismatch_number;
    }

    //!\brief Get mismatch_number.
    constexpr size_type get_mismatch_number() const noexcept
    {
        return _mismatch_number;
    }

    //!\brief Set melting temperature range.
    bool set_Tm_range(float_type min_Tm, float_type max_Tm)
    {
        _primer_melt_range = Range<float_type>{min_Tm, max_Tm};
    }

    //!\brief Set melting temperature range.
    std::pair<float_type, float_type> get_Tm_range()
    {
        return std::make_pair<float_type, float_type>{_primer_melt_range.min, _primer_melt_range.max};
    }

    //!\brief Get minimal melting temperature.
    constexpr float_type get_min_Tm() const noexcept
    {
        return _primer_melt_range.min;
    }

    //!\brief Get maximal melting temperature.
    constexpr float_type get_min_Tm() const noexcept
    {
        return _primer_melt_range.max;
    }

    //!\brief Set molar Natrium concentration if salt adjusted method for primer Tm.
    bool set_Na(float_type Na)
    {
        _Na = Na;
        return true;
    }

    //!\brief Set molar Natrium concentration if salt adjusted method for primer Tm.
    float_type get_Na()
    {
        return Na;
    }

    //!\brief Set method for computing primer melting temperature.
    void set_primer_melt_method(chemistry::method method)
    {
        _primer_melt_method = method;
    }

    //!\brief Get method for computing primer melting temperature.
    chemistry::method get_primer_melt_method()
    {
        return _primer_melt_method;
    }

    //!\brief Set distance_range.
    constexpr bool set_transcript_range(std::pair<size_type, size_type> transcript_range)
    {
        assert(transcript_range.first >= 50 && transcript_range.second <= 2000);
        _distance_range = distance_range;
    }

    //!\brief Get distance_range.
    constexpr std::pair<size_type, size_type> get_transcript_range() const noexcept
    {
        return _transcript_range;
    }

private:
    // lower bound for coverage rate of primer to template match
    float_type _coverage_rate;
    // distance range between forward ending and reverse beginning
    std::pair<size_type, size_type> _distance_range;
    // upper bound for number of mismatches between primer and template
    size_type _mismatch_number;
    //!\brief Constraints for primer length range.
    std::pair<size_type> _transcript_range;
    // root taxonomic identifier
    taxid_type _root_taxid;
    // Natrium concentration
    float_type Na;
};

} // namespace priset
