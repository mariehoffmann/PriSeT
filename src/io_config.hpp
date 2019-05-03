// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================
//          Author: Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
//          Manual: https://github.com/mariehoffmann/PriSeT

#pragma once

#include <cstring>
#include <experimental/filesystem>
#include <iostream>
#include <regex>

#include "errors.hpp"

namespace fs = std::experimental::filesystem;

namespace priset
{

struct io_config
{
private:
    // Directory that contains library, taxonomy, and taxid to accession files.
    fs::path lib_dir;
    // Taxid to accession map in csv format (set by PriSeT).
    fs::path acc_file{};
    // Sequence library file in fasta format (set by PriSeT).
    fs::path fasta_file{};
    // Taxonomy file in csv format (set by PriSeT).
    fs::path tax_file{};
    // Working base directory for storing FM indices, mappings and results.
    fs::path work_dir;
    // Working subdirectory for FM index (set by PriSeT).
    fs::path index_dir;
    // Path to genmap binary (set by PriSeT).
    fs::path genmap_bin;
    // Path to computed FM index (set by PriSeT -> genmap).
    fs::path genmap_idx_dir;
    // Library file extensions.
    std::string ext_fasta = ".fasta";
    std::string ext_acc = ".acc";
    std::string ext_tax = ".tax";

public:
    io_config(fs::path const & lib_dir_, fs::path const & work_dir_) :
        lib_dir{fs::absolute(lib_dir_)},
        work_dir{fs::absolute(work_dir_)},
        index_dir{fs::absolute(work_dir_)}
        {
            // parse library directory and assign paths to the .acc, .fasta, and .tax files
            if (!fs::exists(lib_dir))
                std::cout << "ERROR: " << LIB_DIR_ERROR << std::endl, exit(-1);
            for (auto & p : fs::directory_iterator(lib_dir))
            {
                std::cout << p << std::endl;
                if (p.path().extension().compare(ext_acc) == 0)
                    acc_file = p;
                else if (p.path().extension().compare(ext_fasta) == 0)
                    fasta_file = p;
                else if (p.path().extension().compare(ext_tax) == 0)
                    tax_file = p;
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

            // create working directory if not existing after clearing
            char cmd_mkdir[50];
            if (!fs::exists(work_dir))
            {
                sprintf(cmd_mkdir, "mkdir -p %s", work_dir.c_str());
                if (system(cmd_mkdir))
                    std::cout << "ERROR: " << WRK_DIR_ERROR << std::endl, exit(-1);
            }

            // set output directory for FM index, will be created by genmap
            index_dir += fs::path("/index");
            if (fs::exists(index_dir))
            {
                char cmd_rm[50];
                sprintf(cmd_rm, "rm -r %s", index_dir.c_str());
                if (system(cmd_rm))
                    std::cout << "ERROR: Could not remove index directory " << index_dir << std::endl, exit(0);
            }
        };

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

    // Return path of genmap binary
    fs::path get_genmap_binary() const noexcept
    {
        return genmap_bin;
    }

    // Return directory where FM index is stored
    fs::path get_genmap_idx_dir() const noexcept
    {
        return genmap_idx_dir;
    }

    // Return taxonomy file with absolute path as filesystem::path object.
    fs::path get_tax_file() const noexcept
    {
        return tax_file;
    }

};

} // namespace priset
