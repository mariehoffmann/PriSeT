// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================
//          Author: Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
//          Manual: https://github.com/mariehoffmann/PriSeT

#pragma once

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

#include <type_traits>

#define SEQAN_APP_VERSION "1.0.0"

//#include <seqan/index.h>

// ignore unused variable warnings from https://github.com/cpockrandt/genmap
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"

#include "../submodules/genmap/src/common.hpp"
#include "../submodules/genmap/src/genmap_helper.hpp"
#include "../submodules/genmap/src/indexing.hpp"
#include "../submodules/genmap/src/mappability.hpp"

#include "../submodules/genmap/include/lambda/src/mkindex_saca.hpp"
#include "../submodules/genmap/include/lambda/src/mkindex_misc.hpp"
#include "../submodules/genmap/include/lambda/src/mkindex_algo.hpp"

#pragma GCC diagnostic pop

#include "io_config.hpp"
#include "primer_config.hpp"

namespace priset
{

/*
 * Create FM index with `genmap` binary and store in io_cfg.genmap_idx_dir
 */
template<typename io_config>
int fm_index(io_config & io_cfg)
{
    pid_t pid;
    if ((pid = fork()) == -1)
        std::cout << "ERROR: " << FORK_ERROR << std::endl, exit(0);
    if (pid == 0) {
        execl(io_cfg.get_genmap_binary().c_str(), "genmap", "index", "-F", &io_cfg.get_fasta_file().string()[0u],
            "-I", &io_cfg.get_index_dir().string()[0u], NULL);
        std::cout << "ERROR: " << EXECV_ERROR << std::endl, exit(0);
    }
    else
    {
        wait(NULL);
    }
    return 0;
}

/*
 * Map frequent k-mers to exisiting FM index with `genmap` without file IO
 */
template<typename io_config, typename primer_config>
int fm_map(io_config & io_cfg, primer_config & primer_cfg)
{
    // omit file I/O
    // set flags: genmap map -I io_cfg.genmap_idx_dir -O cfg.genmap_map_dir
    // -K std::string(primer_cfg.get_primer_length_range.min) -E 1 -r -d -fl

    // ./bin/genmap map -I ~/tmp/genmap -O ~/tmp/genmap -K 8 -E 1 --raw -t -d -fl
    std::string argv_str = "-I " + std::string(io_cfg.get_index_dir()) +
                            "-O " + std::string(io_cfg.get_mapping_dir()) +
                            "-K " + std::to_string(primer_cfg.get_primer_length_range().min) + "-E 1 -r -d -fl";
    int argc = 13;
    const char *argv = argv_str.c_str();
    mappabilityMain(argc, &argv);

    //io_cfg.genmap_map_dir += "/genmap_map";

    /*
    // ./bin/genmap map -I ~/tmp/genmap -O ~/tmp/genmap -K 8 -E 1 --raw -t -d -fl
    if ((pid = fork()) == -1)
        std::cout << "ERROR: " << FORK_ERROR << std::endl, exit(0);
    if (pid == 0) {
        std::cout << "start mapping ...\n";
        // create mapping output directory, which has to be already existent for genmap!
        system(&("mkdir " + io_cfg.genmap_map_dir)[0u]);
        // create mapping of k-mers

        execl(io_cfg.genmap_bin.c_str(), "genmap", "map",
            "-I", &io_cfg.genmap_idx_dir[0u],
            "-O", &(cfg.genmap_map_dir)[0u],
            "-K", std::string(primer_cfg.get_primer_length_range.min),
            "-E", "1", "-r", "-d", "-fl", NULL);
        std::cout << "ERROR: " << EXECV_ERROR << std::endl, exit(0);
    }
    else
        wait(NULL);
    */
    return 0;
}

} // namespace priset
