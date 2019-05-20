// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================
//          Author: Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
//          Manual: https://github.com/mariehoffmann/PriSeT

// Interval type for storing upper and lower bounds.

#pragma once

#include <stdlib.h>

#include "primer_config.hpp"

namespace priset
{

/*
 * Datatype to store a kmer as an alphabet sequence, a melting temperature and
 * a unique integer ID.
 */
template<typename primer_config_type>
struct kmer
{
private:
    // Pointer to primer configurator.
    typename std::add_pointer_t<primer_config_type const> primer_cfg{nullptr};
    //primer_config_type &primer_cfg{};
    // Pointer to DNA sequence of k-mer
    //typename std::add_pointer_t<typename primer_config_type::sequence_type const> sequence{nullptr};
    typename primer_config_type::sequence_type const sequence;

    // Unique numerical identifier.
    typename primer_config_type::kmer_ID_type ID;

    //typename primer_config_type::sequence_type sequence;

    // Melting temperature of k-mer sequence.
    typename primer_config_type::float_type Tm;

public:
    constexpr kmer() = default;
    kmer(primer_config_type const & primer_cfg_, typename primer_config_type::sequence_type const & sequence_) :
        primer_cfg(&primer_cfg_), sequence(sequence_)
    {
        Tm = primer_cfg.compute_primer_Tm(sequence_);
    }
    ~kmer() = default;
};

/*
 * Datatype to store matches of two k-mers within a list of accessions.
 */
template<typename kmer_type>
struct match
{
private:

    // k-mer identifier for forward (5') primer sequence
    kmer_type kmer_fwd;
    // k-mer identifier for reverse (3') primer sequence
    kmer_type kmer_rev;
    // taxid of all accessions in the accession list
    typename kmer_type::primer_config_type::taxid_type taxid;
    // absolute difference of their melting temperatures
    typename kmer_type::primer_config_type::float_type Tm_delta;
    // accessions from library where both k-mers match and which are directly assigned
    // to the taxid (e.g. taxid is not the lca)
    // TODO: decide to store here also the location asscociated with an accession
    std::vector<typename kmer_type::primer_cfg::accession_type> accession_list;
public:
    //using = ;
    constexpr match() = default;
    match(kmer_type kmer_fwd, kmer_type kmer_rev)
    {
        Tm_delta = std::abs(kmer_fwd.get_Tm() - kmer_fwd.getTm());
    }
    ~match() = default;
};

// TODO: use key value tuple and use ordered set or map
/*
* Struct to associate a taxonomic identifier with a set of accession numbers and sequences.
*/
template<typename taxid_type=unsigned int, typename accession_type=std::string,  typename sequence_type=std::string>
struct Reference
{
    //!\brief Taxonomic ID associated to the reference collection.
    const taxid_type taxid;

    //!\brief Default constructor.
    constexpr Reference() : taxid(static_cast<taxid_type>(1)) {}
    //!\brief Copy constructor.
    constexpr Reference(Reference const &) = default;
    //!\brief Copy construction via assignment.
    constexpr Reference & operator=(Reference const &) = default;
    //!\brief Move constructor.
    constexpr Reference (Reference &&) = default;
    //!\brief Move assignment.
    constexpr Reference & operator=(Reference &&) = default;
    //!\brief Use default deconstructor.
    ~Reference() = default;

    //!\brief Construct by taxid.
    explicit constexpr Reference(taxid_type const _taxid) noexcept : taxid(_taxid) {}

    void insert(accession_type const accession, sequence_type const sequence)
    {
        assert(accession.size() > 0 && sequence.size() > 0);
        _accessions.push_back(accession);
        _sequences.push_back(sequence);
    }

    //!\brief Erase a sequence given its accession.
    bool erase(accession_type accession)
    {
        auto it = std::find(_accessions.begin(), _accessions.end(), accession);
        if (it != _accessions.end())
        {
            size_t offset = it - _accessions.begin();
            _accessions.erase(it);
            _sequences.erase(_sequences.begin() + offset);
        }
        else
            return false;
        return true;
    }

private:
    //!\brief Accessions associated to taxid.
    std::vector<accession_type> _accessions;
    //!\brief List of sequences associated to the accessions.
    std::vector<sequence_type> _sequences;
};

} // priset
