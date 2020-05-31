// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================
//          Author: Marie Hoffmann <ozymandiaz147 AT gmail.com>
//          Manual: https://github.com/mariehoffmann/PriSeT

// Primer Settings.

#pragma once

#include <iostream>

namespace priset
{

#define WORD_SIZE ((uint8_t)64)

// The number of tailing bits reserved in a KmerID for holding the integer compression of a kmer DNA sequence.
#define CODE_SIZE ((uint8_t)52)

// The number of leading bits reserved in a KmerID to store kmer lengths (= primer_max_length - primer_min_length + 1).
// Maximal possible length is 64 - CODE_SIZE.
#define PREFIX_SIZE ((uint8_t)10)

// Mask selection, i.e. 10 leading bits and rest 0 or ~(1 << 54) - 1)
#define PREFIX_SELECTOR ((uint64_t)(0b1111111111ULL << 54))

// Primer set size we seek to optimize towards coverage.
#define PRIMER_SET_SIZE ((uint8_t)1)

// The minimal possible primer length (or a kmer).
#define PRIMER_MIN_LEN ((uint8_t)16)

// The maximal possible primer length (or a kmer).
#define PRIMER_MAX_LEN ((uint8_t)25)

// The default minimal transcript length.
#define TRANSCRIPT_MIN_LEN ((uint32_t)30)

// The default maximal transcript length.
#define TRANSCRIPT_MAX_LEN ((uint32_t)800)

// The default minimal primer melting temperature. Recommended are 52 °C.
#define PRIMER_MIN_TM ((uint8_t)52)

// The default maximal primer melting temperature. Recommended are 58 °C.
#define PRIMER_MAX_TM ((uint8_t)58)

// Maximal permitted temperature difference [Kelvin] of primers.
// Differences above 5 Kelvin can lead to no amplification and is therefore default.
#define PRIMER_DTM ((uint8_t)5)

// The lower bound for relative CG content. Default is 40 %.
#define CG_MIN_CONTENT ((float).4)

// The upper bound for relative CG content. Default is 60 %.
#define CG_MAX_CONTENT ((float).6)

// The minimal distance (bp) between two identical kmers on same reference.
#define TRAP_DIST ((uint32_t)400)

/* The lower kmer frequency cutoff in percentage, i.e. all kmer occurences
* below will be dropped. Default is 10 %. The absolute value will be computed
* based on the library size.
*/
#define DIGAMMA_PERCENT ((uint8_t)10)

// Lower kmer pair frequency cutoff, i.e. all pair occurences below will be dropped.
#define FREQ_PAIR_MIN ((uint32_t)2)


struct PrimerConfig
{
public:
    // Member types
    using size_type = uint32_t;

    using error_type = uint8_t;

    // interval type for primer length ranges
    using size_interval_type = typename std::pair<size_type, size_type>;

    // The index type to identify kmer sequences.
    using kmer_ID_type = unsigned short int;

private:
    // Number of positions varying from kmer sequence, i.e. number of permitted primer errors.
    error_type E{0};

    // Library size in terms of number of sequences.
    size_type library_size{0};

    // Number of species with at least one reference assigned to it.
    size_type species_count{0};

    uint8_t primer_set_size{PRIMER_SET_SIZE};

    // The minimal primer length. Default is global minimum.
    uint8_t primer_min_len = PRIMER_MIN_LEN;

    // The maximal primer length. Default is global maximum.
    uint8_t primer_max_len = PRIMER_MAX_LEN;

    // The minimal transcript length.
    uint64_t transcript_min_len = TRANSCRIPT_MIN_LEN;

    // The minimal transcript length.
    uint64_t transcript_max_len = TRANSCRIPT_MAX_LEN;

    // The minimal primer melting temperature. Recommended is 52 °C.
    uint8_t primer_min_Tm = PRIMER_MIN_TM;

    // The maximal primer melting temperature.
    uint8_t primer_max_Tm = PRIMER_MAX_TM;

    // Maximal permitted temperature difference [Kelvin] of primers. Recommended is 5 Kelvin.
    uint8_t primer_dTm = PRIMER_DTM;

    // The lower bound for relative CG content. Recommended: 0.4.
    float CG_min_content = CG_MIN_CONTENT;

    // The upper bound for relative CG content. Recommended: 0.6.
    float CG_max_content = CG_MAX_CONTENT;

    // The minimal distance (bp) between two identical kmers on same reference.
    size_type trap_dist = TRAP_DIST;

    // The lower kmer frequency cutoff in percentage.
    float digamma_percent = DIGAMMA_PERCENT;

    /* The lower k-mer frequency cutoff as an absolute value. If set via
    *  set_digamma_absolute, DIGAMMA_PERCENT * |library| will be ignored!
    */
    uint64_t digamma{0};

    // The lower kmer pair frequency cutoff.
    uint64_t digamma_pairs = FREQ_PAIR_MIN;

public:
    // Constructors, destructor and assignment
    // Default constructor.
    PrimerConfig() = default;

    // Default copy constructor.
    constexpr PrimerConfig(PrimerConfig const &) = default;

    // Default copy construction via assignment.
    PrimerConfig & operator=(PrimerConfig const &) = default;

    // Move constructor.
    constexpr PrimerConfig(PrimerConfig && rhs) = default;

    // Move assignment.
    PrimerConfig & operator=(PrimerConfig && rhs) = default;

    // Init with user defined k-mer length and error number
    PrimerConfig(error_type E_)
    {
        E = E_;
    }

    // Use default deconstructor.
    ~PrimerConfig() = default;
    //!\}

    // Set primer set size for which coverage is optimized other than default.
    bool set_primer_set_size(uint8_t const m)
    {
        if (!m)
        {
            std::cerr << "ERROR: primer_set_size needs to be at least one! ";
            std::cerr << "primer_set_size remains " << primer_set_size << std::endl;
            return false;
        }
        primer_set_size = m;
        return true;
    }

    constexpr uint8_t get_primer_set_size() const noexcept
    {
        return primer_set_size;
    }

    /* Set minimal primer length not subceeding PRIMER_MIN_LEN. On success it
    * will be set to l, otherwise remains PRIMER_MIN_LEN.
    */
    bool set_primer_min_len(uint64_t const l)
    {
        if (l < PRIMER_MIN_LEN || l <= PRIMER_MAX_LEN)
        {
            std::cerr << "ERROR: Possible primer length range is [" << PRIMER_MIN_LEN;
            std::cerr << ":" << PRIMER_MAX_LEN << "]. Your request will be ignored!";
            std::cerr << std::endl;
            return false;
        }
        primer_min_len = l;
        if (primer_min_len > primer_max_len)
        {
            std::cerr << "WARNING: primer_min_len > primer_max_len. ";
            std::cerr << "Set primer_max_len to some value larger than primer_min_len!";
            std::cerr << std::endl;
        }
        return true;
    }

    constexpr uint64_t get_primer_min_len() const noexcept
    {
        return primer_min_len;
    }

    /* Set maximal primer length not exceeding PRIMER_MAX_LEN. On success it
    * will be set to l, otherwise remains PRIMER_MAX_LEN.
    */
    bool set_primer_max_len(uint64_t const l)
    {
        if (l > PRIMER_MAX_LEN || l < PRIMER_MIN_LEN)
        {
            std::cerr << "Possible primer length range is [" << PRIMER_MIN_LEN;
            std::cerr << ":" << PRIMER_MAX_LEN << "] °C. Your request will be ignored!";
            std::cerr << std::endl;
            return false;
        }
        if (primer_min_len > primer_max_len)
        {
            std::cerr << "WARNING: primer_min_len > primer_max_len. ";
            std::cerr << "Set primer_min_len to some value smaller than ";
            std::cerr << primer_max_len << "!" << std::endl;
        }
        primer_max_len = l;
        return true;
    }

    constexpr uint64_t get_primer_max_len() const noexcept
    {
        return primer_max_len;
    }

    // Set minimal transcript length other than default.
    bool set_transcript_min_len(size_type const l)
    {
        if (l == 0)
        {
            std::cerr << "Possible transcript length range is [1:";
            std::cerr << ((1UL << 32) - 1) << "]. Your request will be ignored!";
            std::cerr << std::endl;
            return false;
        }
        if (transcript_min_len > transcript_max_len)
        {
            std::cerr << "WARNING: transcript_min_len > transcript_max_len. ";
            std::cerr << "Set transcript_max_len to some value larger than ";
            std::cerr << transcript_max_len << "!" << std::endl;
        }
        transcript_min_len = l;
        return true;
    }

    // Get minimal transcript length.
    constexpr size_type get_transcript_min_len() const noexcept
    {
        return transcript_min_len;
    }

    // Set maximal transcript length other than default.
    bool set_transcript_max_len(size_type const l)
    {
        if (l == 0)
        {
            std::cerr << "ERROR: Possible transcript length range is [1:";
            std::cerr  << ((1UL << 32) - 1) << "]. Your request will be ignored!";
            std::cerr << std::endl;
            return false;
        }
        if (transcript_min_len > transcript_max_len)
        {
            std::cerr << "WARNING: transcript_min_len > transcript_max_len. ";
            std::cerr << "Set transcript_min_len to some value smaller than ";
            std::cerr << transcript_max_len << "!" << std::endl;
        }
        transcript_max_len = l;
        return true;
    }

    // Get maximal transcript length.
    constexpr size_type get_transcript_max_len() const noexcept
    {
        return transcript_max_len;
    }

    // Set minimal primer melting temperature other than default.
    bool set_primer_min_Tm(uint8_t const t)
    {
        if (t < 20 || t > 120)
        {
            std::cerr << "ERROR: Possible melting temperature range is [20:120] °C. ";
            std::cerr << "primer_min_Tm remains " << primer_min_Tm << std::endl;
            return false;
        }
        primer_min_Tm = t;
        if (primer_min_Tm > primer_max_Tm)
        {
            std::cerr << "WARNING: primer_min_Tm > primer_max_Tm. ";
            std::cerr << "Set primer_max_Tm to some value larger than ";
            std::cerr << primer_min_Tm << "!" << std::endl;
        }
        return true;
    }

    // Get minimal primer melting temperature.
    constexpr uint8_t get_primer_min_Tm() const noexcept
    {
        return primer_min_Tm;
    }

    // Set maximal primer melting temperature other than default.
    bool set_primer_max_Tm(uint8_t const t)
    {
        if (t < 20)
        {
            std::cerr << "ERROR: Minimal possible melting temperature is 20 °C. ";
            std::cerr << "primer_min_Tm remains " << primer_min_Tm << std::endl;
            return false;
        }
        primer_min_Tm = t;
        if (primer_min_Tm > primer_max_Tm)
        {
            std::cerr << "WARNING: primer_min_Tm > primer_max_Tm. ";
            std::cerr << "Set primer_max_Tm to some value larger than ";
            std::cerr << primer_min_Tm << "!" << std::endl;
        }
        return true;
    }

    // Get minimal primer melting temperature.
    constexpr uint8_t get_primer_max_Tm() const noexcept
    {
        return primer_max_Tm;
    }

    // Set maximal differences in primer melting temperatures other than default.
    void set_primer_dTm(uint8_t const dTm) noexcept
    {
        primer_dTm = dTm;
    }

    // Get maximal differences in primer melting temperatures.
    constexpr uint8_t get_primer_dTm() const noexcept
    {
        return primer_dTm;
    }

    // Set lower bound for relative CG content other than default.
    void set_CG_min_content(float const CG)
    {
        if (CG > CG_max_content)
        {
            std::cerr << "WARNING: CG_min_content is currently exceeding CG_max_content = ";
            std::cerr << CG_max_content << ". Set CG_max_content to some value larger than ";
            std::cerr << CG << std::endl;
        }
        CG_min_content = CG;
    }

    // Get lower bound for relative CG content other than default.
    constexpr float get_CG_min_content() const noexcept
    {
        return CG_min_content;
    }

    // Set upper bound for relative CG content other than default.
    void set_CG_max_content(float const CG)
    {
        if (CG < CG_min_content)
        {
            std::cerr << "WARNING: CG_max_content is currently subceeding CG_min_content = ";
            std::cerr << CG_min_content << ". Set CG_min_content to some value smaller than ";
            std::cerr << CG << std::endl;
        }
        CG_min_content = CG;
    }

    // Get lower bound for relative CG content other than default.
    constexpr float get_CG_max_content() const noexcept
    {
        return CG_max_content;
    }

    // Set minimal distance (bp) between two identical kmers on same reference other than default.
    void set_trap_dist(size_type const d)
    {
        trap_dist = d;
    }

    // Get minimal distance (bp) between two identical kmers on same reference.
    constexpr size_type get_trap_dist() const noexcept
    {
        return trap_dist;
    }

    // Set size of library in terms of number of sequences/genomes.
    void set_library_size(uint64_t const _library_size)
    {
        library_size = _library_size;
    }

    uint64_t get_library_size() const noexcept
    {
        return library_size;
    }

    // Set number of taxa with references.
    void set_species_count(uint64_t const _species_count)
    {
        species_count = _species_count;
    }

    // Get number of taxa with references.
    uint64_t get_species_count() const noexcept
    {
        return species_count;
    }

    // Set absolute frequency cutoff other than default, which is clade_size*digamma_percent.
    void set_digamma(uint64_t const _digamma)
    {
        digamma = _digamma;
    }

    // Get absolute frequency cutoff for k-mers.
    uint64_t get_digamma() noexcept
    {
        if (!digamma) // if not set, compute based on clade_size
            digamma = size_type(float(digamma_percent)/float(100) * float(species_count));
        return digamma;
    }

    // Set absolut frequency cutoff for k-mer pairs.
    void set_digamma_pairs(uint64_t const _digamma_pairs)
    {
        digamma_pairs = _digamma_pairs;
    }

    // Get absolut frequency cutoff for k-mer pairs.
    constexpr const uint32_t get_digamma_pairs() const noexcept
    {
        return digamma_pairs;
    }

    // Set primer error number (E parameter for genmap's mapping).
    void set_error(error_type const E_)
    {
        assert(E_ < 5);
        E = E_;
    }

    // Get primer error number.
    constexpr error_type get_error() const noexcept
    {
        return E;
    }
};

}  // namespace priset
