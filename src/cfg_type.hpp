#pragma once

#include <string>

#include <seqan/basic.h>


struct cft_type
{

    // The taxonomic identifier type.
    using taxid_type = uint64_t;

    // The type of an accession identifier, e.g. from NCBI GenBank.
    using accession_type = std::string;

};
