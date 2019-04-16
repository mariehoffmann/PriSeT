
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <iostream>
#include <unistd.h>
#include <vector>

#include "errors.hpp"
#include "genmap.hpp"
#include "io_config.hpp"
#include "primer_config.hpp"
#include "taxonomy.hpp"

/*
 * usage        g++ priset.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -o priset
 *              ./priset <path_genmap_bin> <src_dir> <work_dir>
 * e.g.         ./priset $GENMAP ~/tmp/priset/src ~/tmp/priset/work
 *
 * src_dir      path to folder containing fasta (*.fa) and taxonomy file (*.tax)
 * work_dir     path to store indices, mappings, annotations, and other results
 */
int main(int argc, char** argv)
{
    if (argc < 4)
        std::cout << "ERROR: " << ARG_ERROR << std::endl, exit(0);
    // set path prefixes for fasta and taxonomy files
    io_config io_cfg{argv[1], argv[2], argv[3]};
    // build taxonomy in RAM
    taxonomy taxtree{io_cfg.tax_file};

    //set primer constraints and k-mer length
    priset::primer_config<std::string> primer_cfg{};

    run_genmap(io_cfg, primer_cfg);
    std::cout << "fasta_file: " << io_cfg.fasta_file << std::endl;

    // use rotifera.genmap.freq16

    return 0;
}
