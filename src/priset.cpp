// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================
//          Author: Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
//          Manual: https://github.com/mariehoffmann/PriSeT

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <iostream>
#include <unistd.h>
#include <vector>

#include "dna.hpp"
#include "errors.hpp"
#include "fm.hpp"
#include "io_config.hpp"
#include "primer_config.hpp"
#include "taxonomy.hpp"

/*
 * usage        g++ priset.cpp -Wno-write-strings -std=c++17 -lstdc++fs -Wall -Wextra -o priset
 *              ./priset <lib_dir> <work_dir>
 * e.g.         ./src/priset ~/priset/library ~/priset/work
 *
 * src_dir      path to folder containing fasta (*.fa) and taxonomy file (*.tax)
 * work_dir     path to store indices, mappings, annotations, and other results
 */
int main(int argc, char** argv)
{
    if (argc < 3)
        std::cout << "ERROR: " << ARG_ERROR << std::endl, exit(0);

    using sequence_type = std::vector<priset::dna>;

    // set path prefixes for library files
    priset::io_config io_cfg{argv[1], argv[2]};

    // get instance to primer sequence settings
    priset::primer_config<sequence_type> primer_cfg{};

    // build taxonomy in RAM
    priset::taxonomy tax{io_cfg.get_tax_file()};
    tax.print_taxonomy();

    // create FM index
    priset::fm_index(io_cfg);

    // compute k-mer mappings
    priset::fm_map(io_cfg, primer_cfg);

    return 0;
}
