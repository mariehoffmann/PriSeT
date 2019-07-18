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
    std::cout << "MESSAGE: start index recomputation" << std::endl;
    pid_t pid;
    if ((pid = fork()) == -1)
        return FORK_ERROR;
    if (pid == 0) {
        execl(io_cfg.get_genmap_binary().c_str(), "genmap", "index", "-F", &io_cfg.get_fasta_file().string()[0u], "-I", &io_cfg.get_index_dir().string()[0u], NULL);
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
int fm_map(io_cfg_type const & io_cfg, primer_cfg_type const & primer_cfg, TLocations & locations)
{
    std::cout << "STATUS: run genmap::mappability with E = " << primer_cfg.get_error() << std::endl;
    std::cout << "INFO: K = " << primer_cfg.get_primer_length_range().first << std::endl;
    std::string s1 = io_cfg.get_index_dir().string();
    std::string s2 = io_cfg.get_mapping_dir().string();
    // Remark: csv flag triggers `csvComputation` and therefore the population of the (TLocations) locations vector!
    char const * argv[11] = {"map", "-I", s1.c_str(), "-O", s2.c_str(), "-K", std::to_string(primer_cfg.get_primer_length_range().first).c_str(), "-E", std::to_string(primer_cfg.get_error()).c_str(), "--csv", "-fl"};

    mappabilityMain<TLocations>(11, argv, locations);
    return 0;
}

} // namespace priset
