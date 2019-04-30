// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================

/*!\file
* \brief Interval type for storing upper and lower bounds.
* \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
*/

#pragma once

namespace priset
{
template<typename value_type>
struct interval
{
    bool in(value_type const value) const
    {
        assert(min <= max && "max set to a value less than min!");
        return min <= value && value <= max;
    }
    value_type min;
    value_type max;
};

template<typename primer_config_type>
struct kmer
{
private:
    // Unique numerical identifier
    unsigned short int ID;

    // DNA sequence of k-mer
    typename primer_config_type::sequence_type sequence;

    // Its melting temperature
    typename primer_config::float_type Tm;

public:
    kmer(primer_config_type::sequence_type sequence_)
    {
        Tm = io_cfg.compute_primer_Tm();
    }
};

template<typename primer_config_type>
struct match
{
private:
    // k-mer identifier for forward (5') primer sequence
    unsigned short int kmer_ID1;
    // k-mer identifier for reverse (3') primer sequence
    unsigned short int kmer_ID2;
    // taxid of all accessions in the accession list
    typename primer_config_type::taxid_type taxid;
    // difference of their melting temperatures
    typename primer_config_type::float_type Tm_diff;
    // accessions from library where both k-mers match and which are directly assigned
    // to the taxid (e.g. taxid is not the lca)
    std::vector<typename primer_cfg::accession_type> accession_list;
public:
};

} // priset
