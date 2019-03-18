#pragma once

#include <cstring>
#include <iostream>

struct io_config
{
    // path to genmap binary
    std::string genmap_bin;
    // source directory containing fasta and taxonomy file
    std::string src_dir;
    // working directory for storing FM indices and mappings
    std::string work_dir;
    std::string fasta_file;
    std::string tax_file;
    std::string genmap_idx_dir;
    std::string genmap_map_dir;
    char suffix_fasta[4] = ".fa";
    char suffix_tax[5] = ".tax";

    io_config(std::string const genmap_bin_, std::string const src_dir_, std::string const work_dir_) :
        genmap_bin{genmap_bin_}, src_dir{src_dir_}, work_dir{work_dir_}, fasta_file{src_dir_}, tax_file{src_dir_},
        genmap_idx_dir{work_dir_}, genmap_map_dir{work_dir_}
        {};
};
