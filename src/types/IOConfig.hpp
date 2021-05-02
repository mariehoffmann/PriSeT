// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================
//          Author: Marie Hoffmann <ozymandiaz147 AT gmail.com>
//          Manual: https://github.com/mariehoffmann/PriSeT

#pragma once

#include <cstring>
#include <experimental/filesystem>
#include <iostream>
#include <regex>
#include <unordered_set>
#include <stdexcept>

#include "Errors.hpp"
#include "simple_types.hpp"

namespace fs = std::experimental::filesystem;

namespace priset
{

#define MAX_PATH_LENGTH 100

// TODO: rename to io_configurator, instance are then called io_configuration
struct IOConfig
{

public:
    // Default constructor.
    IOConfig() = default;

    // Default copy constructor.
    IOConfig(IOConfig const &) = default;

    // Default copy construction via assignment.
    IOConfig & operator=(IOConfig const &) = default;

    // Move constructor.
    IOConfig(IOConfig && rhs) = default;

    // Move assignment.
    IOConfig & operator=(IOConfig && rhs) = default;

    // Set library and working directory paths.
    // skip_idx = 0  compute FM-index
    // skip_idx = 1  else skip FM-index computation and use index in given directory
    void assign(fs::path const & lib_dir_, fs::path const & work_dir_, bool const _skip_idx = 1)
    {
        lib_dir = fs::canonical(lib_dir_);
        work_dir = fs::canonical(work_dir_);
        if (!fs::exists(work_dir_))
        {
            std::cout << "work dir does not exist, creating ...\n";
            if (!fs::exists(work_dir_.parent_path()))
                fs::create_directory(work_dir_.parent_path());
            fs::create_directory(work_dir_);
        }

        skip_idx = _skip_idx;
        index_dir = work_dir;
        mapping_dir = work_dir;

        if (!fs::exists(lib_dir))
            std::cout << "ERROR: " << LIB_DIR_ERROR << std::endl, exit(-1);
        for (auto & p : fs::directory_iterator(lib_dir))
        {
            if (p.path().extension().compare(ext_acc) == 0)
                acc_file = p;
            else if (p.path().extension().compare(ext_fasta) == 0)
            {
                fasta_file = p;
                char const start = '>';
                std::ifstream stream(fasta_file.string().c_str(), std::ios::in);
                std::string line;
                while (std::getline(stream, line))
                {
                    if (line.at(0) == start)
                        ++library_size;
                }
                std::cout << "INFO\tlibrary size: " << library_size << std::endl;
            }
            else if (p.path().extension().compare(ext_tax) == 0)
                tax_file = p;
        }

        if (!fasta_file.has_filename())
            std::cout << "ERROR: Unable to locate fasta file in: " << lib_dir << std::endl, exit(-1);
        std::cout << "STATUS\tSet fasta file: \t" << fasta_file << std::endl;

        // create working directory if not existing after clearing
        if (!fs::exists(work_dir))
        {
            char cmd[50];
            sprintf(cmd, "mkdir -p %s", work_dir.c_str());
            if (system(cmd))
                std::cout << "ERROR: " << WRK_DIR_ERROR << std::endl, exit(-1);
        }

        // set output directory for FM index, will be created by genmap
        index_dir /= fs::path("index");
        std::cout << "STATUS\tSet index directory: \t" << index_dir << std::endl;

        if (skip_idx && !fs::exists(index_dir))
        {
            std::cerr << "ERROR: You configured PriSeT to skip index computation, but index_dir (" << index_dir << ") does not exists!" << std::endl;
            exit(-1);
        }

        if (!skip_idx && fs::exists(index_dir))
        {
            char cmd[50];
            sprintf(cmd, "rm -r %s", index_dir.c_str());
            if (system(cmd))
                std::cout << "ERROR: Could not remove index directory " << index_dir << std::endl, exit(-1);
        }

        // set output directory for FM index mapping
        mapping_dir /= fs::path("/mapping");
        std::cout << "STATUS\tSet mapping_dir: " << mapping_dir << std::endl;

        // create table output directory
        fs::path result_path = work_dir / "table";
        if (!fs::exists(result_path))
        {
            char cmd[50];
            sprintf(cmd, "mkdir %s", result_path.c_str());
            if (system(cmd))
                std::cout << "ERROR: Creating result table directory = " << result_path << std::endl, exit(-1);
        }
        result_file = result_path / "result.csv";
        primer_info_file = result_path / "primer_info.csv";

        /* Parse accession, id, and taxonomy files to build species set,
        * sequence to accession, and taxid to accessions map. Call first
        * build_species_set to store only taxids of species (and not of higher
        * order).
        */
        build_species_set();
        build_taxid2accs_map();
        build_seqNo2acc_map();
    };

    // Destructor.
    ~IOConfig() = default;

    // Return skip_idx.
    bool get_skip_idx() const noexcept
    {
        return skip_idx;
    }

    // Return accession file with absolute path as filesystem::path object.
    fs::path get_acc_file() const noexcept
    {
        return acc_file;
    }

    // Return library file with absolute path as filesystem::path object.
    fs::path get_fasta_file() const noexcept
    {
        return fasta_file;
    }

    // Return directory where FM index is stored.
    fs::path get_index_dir() const noexcept
    {
        return index_dir;
    }

    // return path to concatenation of text corpus stored in two files index.txt.concat and index.txt.limits
    fs::path get_index_txt_path() const noexcept
    {
        return index_dir / "index.txt";
    }

    // Return directory for FM mappility results.
    fs::path get_mapping_dir() const noexcept
    {
        return mapping_dir;
    }

    // Return primer info file.
    fs::path get_primer_info_file() const noexcept
    {
        return primer_info_file;
    }

    // Return file to store results in csv format.
    fs::path get_result_file() const noexcept
    {
        return result_file;
    }

    // Return taxonomy file with absolute path as filesystem::path object.
    fs::path get_tax_file() const noexcept
    {
        return tax_file;
    }

    // Return working directory.
    fs::path get_work_dir() const noexcept
    {
        return work_dir;
    }

    // Return the library size, i.e. number sequences in FASTA file.
    constexpr uint32_t get_library_size() const noexcept
    {
        return library_size;
    }

    // Return number of species.
    size_t get_species_count() const noexcept
    {
        return species_set.size();
    }

    // Test if given taxonomic ID corresponds to a species.
    bool is_species(Taxid const & taxid) const noexcept
    {
        return (species_set.find(taxid) == species_set.end()) ? 0 : 1;
    }

    // Get next species stored in set. We use the fact that order of unordered_set
    // is fixed § 23.2.5 [unord.req] once no insertion takes place.
    Taxid get_next_species() noexcept
    {
        if (it_species == species_set.cend())
        {
            it_species = species_set.cbegin();
            return Taxid{0};
        }
        Taxid taxid = *it_species;
        ++it_species;
        return taxid;
    }

    // Given a sequence number return its accession ID based on FASTA file.
    Accession get_acc_by_seqNo(TSeqNo const seqNo) const
    {
        if (seqNo2acc_map.count(seqNo))
            return seqNo2acc_map.at(seqNo);
        return "";
    }

    // Given an accession return the taxonomic ID.
    Taxid get_taxid_by_accession(Accession acc) const
    {
        return acc2taxid_map.at(acc);
    }

    // Get taxid given an uncompressed sequence number.
    Taxid get_taxid_by_seqNo(TSeqNo seqNo) const
    {
        assert(seqNo2acc_map.count(seqNo));
        return acc2taxid_map.at(seqNo2acc_map.at(seqNo));
    }

private:
    // Directory that contains library, taxonomy, and taxid to accession files.
    fs::path lib_dir;

    // Working base directory for storing FM indices, mappings and results.
    fs::path work_dir;

    // Flag for indicating if index computation shall be skipped (because it already exists).
    bool skip_idx{1};

    // Taxid to accession map in csv format (set by PriSeT).
    fs::path acc_file{};

    // Sequence library file in fasta format (set by PriSeT).
    fs::path fasta_file{};

    // Taxonomy file in csv format (set by PriSeT).
    fs::path tax_file{};

    // Working subdirectory for FM index (set by PriSeT).
    fs::path index_dir{};

    // Working subdirectory for FM index mappings (set by PriSeT).
    fs::path mapping_dir{};

    // Library file extensions.
    std::string ext_fasta = ".fasta";
    std::string ext_acc = ".acc";
    std::string ext_tax = ".tax";
    std::string ext_id = ".id";

    // Library size in terms of number of accessions (= fasta entries)
    uint64_t library_size{0};

    // Species extracted from tax_file.
    std::set<Taxid> species_set;

    // Iterator for species set.
    std::set<Taxid>::const_iterator it_species = species_set.cbegin();

    // Taxid to accessions map.
    std::unordered_map<Taxid, std::vector<Accession>> taxid2accs_map;

    // For accession/reference ID store assigned taxid.
    std::unordered_map<Accession, Taxid> acc2taxid_map;

    // Association between sequence and accession identifiers.
    std::unordered_map<TSeqNo, Accession> seqNo2acc_map;

    // Path to store result tables to load in Shiny.
    fs::path result_file{};

    // Path to primer info file (sequences and chemical attributes)
    fs::path primer_info_file{};

    // Fill species set based on taxonomy file with row format taxid,p_taxid,is_species.
    void build_species_set()
    {
        species_set.clear();
        char const delim = ',';
        std::ifstream stream(tax_file.string().c_str(), std::ios::in);
        std::string row;
        // count listed species in tax file with row format "taxid,parent_taxid,is_species"
        std::string const is_species = "1";
        std::getline(stream, row); // ignore header line
        while (std::getline(stream, row))
        {
            if (!row.compare(row.size() - 1, 1, is_species))
            {
                Taxid taxid = std::stoi(row.substr(0, row.find(delim)));
                species_set.insert(taxid);
            }
        }
        stream.close();
        // update set iterator
        it_species = species_set.cbegin();
    }

    // Fill taxid to accessions map
    void build_taxid2accs_map()
    {
        acc2taxid_map.clear();
        taxid2accs_map.clear();
        char const delim = ',';
        std::ifstream stream(acc_file.string().c_str(), std::ios::in);
        std::string row;
        std::getline(stream, row); // ignore header
        while (std::getline(stream, row))
        {
            size_t pos1{0};
            size_t pos2 = row.find(delim, pos1 + 1);
            Taxid taxid = std::stoi(row.substr(pos1, pos2));
            if (species_set.find(taxid) == species_set.end())
                continue;
            std::vector<Accession> accs;
            pos1 = pos2;
            pos2 = row.find(delim, pos1 + 1);
            while (pos1 != std::string::npos)
            {
                Accession acc = row.substr(pos1 + 1, pos2 - pos1 - 1);
                accs.push_back(acc);
                pos1 = pos2;
                pos2 = row.find(delim, pos1 + 2);
                acc2taxid_map[acc] = taxid;
            }
            taxid2accs_map[taxid] = accs;
        }
        stream.close();
    }

    // Helper routine to build seqNo2acc_map.
    // Extract accession ID from FASTA file and use previously built acc2taxid_map
    // to determine taxid. If taxid corresponds to a species it is written into
    // dictionary.
    bool build_seqNo2acc_map()
    {
        seqNo2acc_map.clear();
        int seqNo{0};  // 1-based sequence number
        std::unordered_set<std::string> accs_seen;
        char const start = '>';
        char const delim = ' ';
        std::ifstream stream(fasta_file.string().c_str(), std::ios::in);
        std::string row;
        size_t pos;
        while (std::getline(stream, row))
        {
            if (row[0] != start)
                continue;
            ++seqNo;
            pos = row.find(delim);
            Accession acc = row.substr(1, pos-1);
            if (accs_seen.find(acc) != accs_seen.end())
                throw std::domain_error("ERROR: duplicate accession ID (" + acc + ") fasta file");
            accs_seen.insert(acc);
            try
            {
                Taxid taxid = acc2taxid_map[acc];
                if (species_set.find(taxid) != species_set.end())
                {
                    seqNo2acc_map[seqNo] = acc;
                }
            }
            catch (const std::domain_error & e)
            {
                stream.close();
                std::cerr << "ERROR: extracted accession from FASTA file is not listed in related accession file!" << std::endl;
            }
        }
        stream.close();
        return (seqNo2acc_map.size()) ? true : false;
    }

};

// The NULL value of an I/O Configurator.
IOConfig NULL_IOConfig = IOConfig{};

} // namespace priset
