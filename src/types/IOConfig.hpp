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

// #include "chemistry.hpp"
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
    void assign(fs::path const & lib_dir_, fs::path const & work_dir_, bool const idx_only_flag_, bool const skip_idx_flag_)
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

        idx_only_flag = idx_only_flag_;
        skip_idx_flag = skip_idx_flag_;
        index_dir = work_dir;
        mapping_dir = work_dir;
        // TODO: path to PriSeT git repos as argument
        // genmap_bin = "~/git/PriSet_git2/PriSeT/submodules/genmap/bin/genmap";

        if (!fs::exists(lib_dir))
            std::cout << "ERROR: " << LIB_DIR_ERROR << std::endl, exit(-1);
        for (auto & p : fs::directory_iterator(lib_dir))
        {
            //std::cout << p << std::endl;
            if (p.path().extension().compare(ext_acc) == 0)
                acc_file = p;
            else if (p.path().extension().compare(ext_fasta) == 0)
                fasta_file = p;
            else if (p.path().extension().compare(ext_tax) == 0)
                tax_file = p;
            else if (p.path().extension().compare(ext_id) == 0)
            {
                id_file = p;
                // set library size
                // faster with  wc -l root_10190.id
                //std::system("ls -l >test.txt"); // execute the UNIX command "ls -l >test.txt"
                //std::cout << std::ifstream("test.txt").rdbuf();

                FILE * infile = fopen(id_file.string().c_str(), "r");
                int c;
                while (EOF != (c = getc(infile)))
                    if (c == '\n')
                        ++library_size;
                std::cout << "INFO: library size = " << library_size << std::endl;
            }
        }
        if (!acc_file.has_filename())
            std::cout << "ERROR: Unable to locate accession file in: " << lib_dir << std::endl, exit(-1);


        std::cout << "STATUS\tSet accessions file: \t" << acc_file << std::endl;
        if (!fasta_file.has_filename())
            std::cout << "ERROR: Unable to locate fasta file in: " << lib_dir << std::endl, exit(-1);

        std::cout << "STATUS\tSet fasta file: \t" << fasta_file << std::endl;
        if (!tax_file.has_filename())
            std::cout << "ERROR: Unable to locate taxonomy file in: " << lib_dir << std::endl, exit(-1);
        std::cout << "STATUS\tSet taxonomy file: \t" << tax_file << std::endl;
        if (!id_file.has_filename())
            std::cout << "ERROR: Unable to locate id file in: " << lib_dir << std::endl, exit(-1);
        std::cout << "STATUS\tSet id file: \t" << id_file << std::endl;

        // create working directory if not existing after clearing
        if (!fs::exists(work_dir))
        {
            char cmd[50];
            sprintf(cmd, "mkdir -p %s", work_dir.c_str());
            if (system(cmd))
                std::cout << "ERROR: " << WRK_DIR_ERROR << std::endl, exit(-1);
        }

        // set output directory for FM index, will be created by genmap
        index_dir /= fs::path("/index");
        if (skip_idx_flag && !fs::exists(index_dir))
        {
            std::cerr << "ERROR: Index computation flag is set to 0, but index_dir (" << index_dir << ") does not exists!" << std::endl;
            exit(-1);
        }

        if (!skip_idx_flag && fs::exists(index_dir))
        {
            char cmd[50];
            sprintf(cmd, "rm -r %s", index_dir.c_str());
            if (system(cmd))
                std::cout << "ERROR: Could not remove index directory " << index_dir << std::endl, exit(0);
        }
        // set output directory for FM index mapping
        mapping_dir /= fs::path("/mapping");
        std::cout << "set mapping_dir = " << mapping_dir << std::endl;

        // create table output directory
        fs::path result_path = work_dir / "table";
        if (!fs::exists(result_path))
        {
            char cmd[50];
            sprintf(cmd, "mkdir %s", result_path.c_str());
            if (system(cmd))
                std::cout << "ERROR: Creating result table directory = " << result_path << std::endl, exit(-1);
        }
        result_file = result_path / "results.csv";
        primer_info_file = result_path / "primer_info.csv";
        script_file = get_work_dir() / "app" / "app.R";
        std::cout << "STATUS\tSet R script file: " << script_file << std::endl;
        script_runner = get_work_dir() / "app" / "app_run.R";

        if (!fs::exists(script_runner))
        {
            std::ofstream ofs;
            ofs.open(script_runner);
            ofs << "library(shiny)\nrunApp(" << script_file << ")\n";
            ofs.close();
        }

        /* Parse accession and taxonomy file to build species set and taxid to
        * accessions map. Call first build_species_set to store only taxids of
        * species (and not of higher order).
        */
        build_species_set();
        build_taxid2accs_map();
    };

    // Destructor.
    ~IOConfig() = default;

    // Return skip_idx flag.
    bool idx_only() const noexcept
    {
        return idx_only_flag;
    }

    // Return skip_idx flag.
    bool skip_idx() const noexcept
    {
        return skip_idx_flag;
    }

    // Return accession file with absolute path as filesystem::path object.
    fs::path get_acc_file() const noexcept
    {
        return acc_file;
    }

    fs::path get_id_file() const noexcept
    {
        return id_file;
    }

    constexpr uint32_t get_library_size() const noexcept
    {
        return library_size;
    }

    // Return number of species.
    size_t get_species_count() const noexcept
    {
        return species_set.size();
    }

    // Return template file for shiny app.
    fs::path get_app_template() const noexcept
    {
        return app_template;
    }

    // Return library file with absolute path as filesystem::path object.
    fs::path get_fasta_file() const noexcept
    {
        return fasta_file;
    }

    // Return path of genmap binary
    fs::path get_genmap_binary() const noexcept
    {
        return work_dir / "genmap_bin";
    }

    // Return directory where FM index is stored
    fs::path get_index_dir() const noexcept
    {
        return index_dir;
    }

    // path to index plus basename without suffix, i.e. <path_to_index>/<basename>
    fs::path get_index_base_path() const noexcept
    {
        return index_dir / "index";
    }

    // path to index plus basename without suffix, i.e. <path_to_index>/<basename>
    fs::path get_index_base_path_ids() const noexcept
    {
        return index_dir / "index.ids";
    }

    // return path to concatenation of text corpus stored in two files index.txt.concat and index.txt.limits
    fs::path get_index_txt_path() const noexcept
    {
        return index_dir / "index.txt";
    }

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

    // Return path to R script file to be run in terminal.
    fs::path get_script_file() const noexcept
    {
        std::cout << "Enter get_script_file, " << script_file << std::endl;
        return script_file;
    }

    // Return path to R script starter file.
    fs::path get_script_runner() const noexcept
    {
        return script_runner;
    }

    // Return taxonomy file with absolute path as filesystem::path object.
    fs::path get_tax_file() const noexcept
    {
        return tax_file;
    }

    fs::path get_work_dir() const noexcept
    {
        return work_dir;
    }

private:
    // Directory that contains library, taxonomy, and taxid to accession files.
    fs::path lib_dir;
    // Working base directory for storing FM indices, mappings and results.
    fs::path work_dir;
    // Flag for indicating if index computation shall be skipped (because it already exists).
    bool skip_idx_flag{0};
    // Flag for indicating to do index computation exclusively.
    bool idx_only_flag{0};
    // Taxid to accession map in csv format (set by PriSeT).
    fs::path acc_file{};
    // Sequence library file in fasta format (set by PriSeT).
    fs::path fasta_file{};

    // 1-based IDs and accession numbers.
    fs::path id_file{};
    // Taxonomy file in csv format (set by PriSeT).
    fs::path tax_file{};
    // Working subdirectory for FM index (set by PriSeT).
    fs::path index_dir;
    // Working subdirectory for FM index mappings (set by PriSeT).
    fs::path mapping_dir;

    // Library file extensions.
    std::string ext_fasta = ".fasta";
    std::string ext_acc = ".acc";
    std::string ext_tax = ".tax";
    std::string ext_id = ".id";

    // Library size in terms of number of accessions (= fasta entries)
    uint64_t library_size{0};

    // Species extracted from tax_file.
    std::unordered_set<Taxid> species_set;

    // Taxid to accessions map.
    std::unordered_map<Taxid, std::vector<Accession>> taxid2accs_map;

    // For accession/reference ID store assigned taxid.
    std::unordered_map<Accession, Taxid> acc2taxid_map;

    // Path to R shiny app template
    fs::path app_template = "../PriSeT/src/app_template.R";

    // R script for launching shiny app.
    fs::path script_runner;

    // Path to generated copy of R script to be run in terminal.
    fs::path script_file;

    // Path to store result tables to load in Shiny.
    fs::path result_file;

    // Path to primer info file (sequences and chemical attributes)
    fs::path primer_info_file;

    // Fill species set based on taxonomy file with row format taxid,p_taxid,is_species.
    void build_species_set()
    {
        char const delim = ',';
        std::ifstream stream(tax_file.string().c_str(), std::ios::in);
        std::string row;
        // count listed species in tax file with row format "taxid,parent_taxid,is_species"
        std::string const is_species = "1";
        std::getline(stream, row); // ignore header line
        while (std::getline(stream, row))
        {
            if (row.compare(row.size() - 2, 1, is_species) == 1)
            {
                Taxid taxid = std::stoi(row.substr(0, row.find(delim)));
                species_set.insert(taxid);
            }
        }
    }

    // Fill taxid to accessions map
    // TODO: write test
    void build_taxid2accs_map()
    {
        char const delim = ',';
        std::ifstream stream(acc_file.string().c_str(), std::ios::in);
        std::string row;
        std::getline(stream, row); // ignore header
        while (std::getline(stream, row))
        {
            size_t pos1{0};
            size_t pos2 = row.find(delim, pos1);
            Taxid taxid = std::stoi(row.substr(pos1, pos2));
            if (species_set.find(taxid) == species_set.end())
                continue;
            std::vector<Accession> accs;
            pos1 = pos2;
            pos2 = row.find(delim, pos1);
            while (pos1 != std::string::npos)
            {
                Accession acc = row.substr(pos1, pos2);
                accs.push_back(acc);
                pos1 = pos2;
                pos2 = row.find(delim, pos1);
                acc2taxid_map[acc] = taxid;
            }
            taxid2accs_map[taxid] = accs;
        }
    }

};

} // namespace priset