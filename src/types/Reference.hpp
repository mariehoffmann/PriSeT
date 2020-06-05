// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================

/*!\file
* \brief Reference type for a taxid and its reference sequences.
* \author Marie Hoffmann <ozymandiaz147 AT gmail.com>
*/

#pragma once

#include <algorithm>
#include <string>
#include <vector>

namespace priset
{
// TODO: use key value tuple and use ordered set or map
/*
* Struct to associate a taxonomic identifier with a set of accession numbers and sequences.
*/
template<typename taxid_type=unsigned int, typename Accession=std::string,  typename sequence_type=std::string>
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

    void insert(Accession const accession, sequence_type const sequence) noexcept
    {
        _accessions.push_back(accession);
        _sequences.push_back(sequence);
    }

    //!\brief Erase a sequence given its accession.
    bool erase(Accession accession)
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
    std::vector<Accession> _accessions;
    //!\brief List of sequences associated to the accessions.
    std::vector<sequence_type> _sequences;
};

} // priset
