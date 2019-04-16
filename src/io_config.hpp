#pragma once

#include <cstring>
#include <filesystem>
#include <iostream>
#include <regex>

struct io_config
{
private:
    // path to genmap binary
    std::string genmap_bin;
    // source directory containing fasta and taxonomy file
    std::string src_dir;
    // working directory for storing FM indices and mappings
    std::filesystem::path work_dir;
    std::string fasta_file;
    std::string tax_file;
    std::filesystem::path genmap_idx_dir;
    std::filesystem::path genmap_map_dir;
    char suffix_fasta[4] = ".fa";

    std::regex suffix_fasta_rx("\\.(fa|fasta)");

    char suffix_tax[5] = ".tax";

public:
    io_config(std::string const genmap_bin_, std::filesystem::path const & src_dir_, std::filesystem::path const & work_dir_) :
        genmap_bin{genmap_bin_}, src_dir{std::filesystem::absolute(src_dir_)},
        work_dir{std::filesystem::absolute(work_dir_)},
        genmap_idx_dir{std::filesystem::absolute(work_dir_)},
        genmap_map_dir{std::filesystem::absolute(work_dir_)}
        {
            // todo: locate fasta and taxonomy files here, move code from genmap.hpp
            // check if source directory exists and contains one *.fa and one *.tax file
            DIR *dirp;
            struct dirent *dp;
            std::cout << io_cfg.fasta_file << std::endl;
            if ((dirp = opendir(io_cfg.fasta_file.c_str())) == NULL)
                std::cout << "ERROR: " << SRC_DIR_ERROR << std::endl, exit(0);

            bool fasta_set = false, tax_set = false;
            for (const auto & entry : fs::directory_iterator(io_cfg.get_fasta_file()))
            {
                std::cout << entry.path() << std::endl;
                if (!fasta_set && strstr(entry, io_cfg.get_suffix_fasta())
                    io_cfg.set_fasta_file(entry);
                else if (!tax_set && strstr(entry, io_cfg.get_suffix_tax())
                    io_cfg.set_tax_file(entry);
            }

            do
            {
                errno = 0;
                if ((dp = readdir(dirp)) != NULL)
                {
                    // TODO: use regex to match common fasta file suffixes
                    if (!fasta_set && dp->d_reclen > 3 && std::strstr(dp->d_name, io_cfg.suffix_fasta))
                        io_cfg.fasta_file += "/" + std::string(dp->d_name), fasta_set = true;
                    else if (!tax_set && dp->d_reclen > 4 && strstr(dp->d_name, io_cfg.suffix_tax))
                        io_cfg.tax_file += "/" + std::string(dp->d_name), tax_set = true;
                }
                else
                {
                    if (errno == 0) {
                        closedir(dirp);
                        std::cout << "ERROR: " << NO_SRC_ERROR << std::endl, exit(0);
                    }
                    closedir(dirp);
                    std::cout << "ERROR: " << SRC_READ_ERROR << std::endl, exit(0);
                }
            }
            while (dp != NULL && !fasta_set && !tax_set);

            std::cout << "Found fasta file:\t" << io_cfg.fasta_file << "\nFound taxonomy file:\t" << io_cfg.tax_file << std::endl;
        };

    // setter and getter for fasta file
    void set_fasta_file(std::string fasta_file_)
    {
        fasta_file{fasta_file_};
    }

    std::path get_fasta_file()
    {
        return src_dir += std::filesystem::path(fasta_file);
    }

    void set_tax_file(std::string tax_file_)
    {
        tax_file{tax_file_};
    }

    std::path get_tax_file()
    {
        return src_dir += std::filesystem::path(tax_file);
    }


    // return fasta file with absolute path as filesystem::path object
    std::path get_tax_file()
    {
        return src_dir += std::filesystem::path(tax_file);
    }

};
