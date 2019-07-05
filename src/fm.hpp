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

// TODO: place proper seqan version
#define SEQAN_APP_VERSION "1.0.0"

//#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>

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

#include "io_cfg_type.hpp"
#include "primer_cfg_type.hpp"
#include "types.hpp"
#include "utilities.hpp"

namespace priset
{
/*
 * Create FM index with `genmap` binary and store in io_cfg.genmap_idx_dir
 */
int fm_index(io_cfg_type const & io_cfg)
{
    if (io_cfg.get_skip_idx())
        return;
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
 * io_cfg_type              I/O configurator type
 * primer_cfg_type          primer configurator type
 * TLocations               type for storing locations
 * TDirectoryInformation    directory information type
 */
int fm_map(io_cfg_type const & io_cfg, primer_cfg_type const & primer_cfg, TLocations & locations, TDirectoryInformation & directoryInformation)
{
    std::cout << "/Users/troja/priset/335928/work/index == ? " << io_cfg.get_index_dir() << std::endl;
    std::cout << "/Users/troja/priset/335928/work/mapping == ? " << io_cfg.get_mapping_dir() << std::endl;
    std::string s1 = io_cfg.get_index_dir().string();
    std::string s2 = io_cfg.get_mapping_dir().string();
    std::cout << "s1 = " << s1 << std::endl;
    std::cout << "s2 = " << s2 << std::endl;

    //char const * argv[12] = {"map", "-I", io_cfg.get_index_dir().c_str(), "-O", io_cfg.get_mapping_dir().c_str(), "-K", "18", "-E", "0", "--raw", "-fl", NULL};
//    char const * argv[12] = {"map", "-I", s1.c_str(), "-O", s2.c_str(), "-K", "8", "-E", std::to_string(primer_cfg.get_error()).c_str(), "--raw", "-fl", NULL};
    char const * argv[12] = {"map", "-I", s1.c_str(), "-O", s2.c_str(), "-K", "7", "-E", std::to_string(2).c_str(), "--raw", "-fl", NULL};


    /*
    argv[0] = map
    argv[1] = -I
    argv[2] = /Users/troja/priset/335928/work/index
    argv[3] = -O
    argv[4] = /Users/troja/priset/335928/work/mapping
    argv[5] = -K
    argv[6] = 18
    argv[7] = -E
    argv[8] = 0
    argv[9] = --raw
    argv[10] = -fl
        */
    mappabilityMain(11, argv);
    return 0;
}

} // namespace priset
