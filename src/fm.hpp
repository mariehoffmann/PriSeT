// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================
//          Author: Marie Hoffmann <ozymandiaz147 AT gmail.com>
//          Manual: https://github.com/mariehoffmann/PriSeT

#pragma once

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <sys/types.h>
#include <sys/wait.h>
#include <type_traits>
#include <unistd.h>
#include <vector>

// TODO: place proper seqan version
#define SEQAN_APP_VERSION "1.0.0"

//#include <seqan/arg_parse.h>
#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/index.h>

// ignore unused variable warnings from https://github.com/cpockrandt/genmap
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"

#define GENMAP_DEBUG 0

#include "../submodules/genmap/include/lambda/src/mkindex_saca.hpp"
#include "../submodules/genmap/include/lambda/src/mkindex_misc.hpp"
#include "../submodules/genmap/include/lambda/src/mkindex_algo.hpp"
#include "../submodules/genmap/src/common.hpp"
#include "../submodules/genmap/src/genmap_helper.hpp"
#include "../submodules/genmap/src/indexing.hpp"
#include "../submodules/genmap/src/mappability.hpp"

#pragma GCC diagnostic pop

#include "types/all.hpp"
#include "utilities.hpp"

namespace priset
{
/*
 * Create FM index with `genmap` binary and store in io_cfg.genmap_idx_dir
 */
// int fm_index2(IOConfig const & io_cfg)
// {
//     std::cout << "MESSAGE: start index recomputation" << std::endl;
//     std::string cmd = io_cfg.get_genmap_binary().string() + " index -F " + io_cfg.get_fasta_file().string() + " -I " + io_cfg.get_index_dir().string();
//     std::system(cmd.c_str());
//     return 0;
// }

int fm_index(IOConfig const & io_cfg)
{
    int const argc = 5;
    std::string const fasta_file = io_cfg.get_fasta_file().string();
    char const * argv[argc] = {"index",
        "-F", fasta_file.c_str(),
        "-I", io_cfg.get_index_dir().c_str()};
    indexMain(argc, argv);
    return 0;
}

/*
 * Map frequent k-mers to exisiting FM index with `genmap` without file IO
 * IOConfig              I/O configurator type
 * PrimerConfig          primer configurator type
 * TLocations               type for storing locations
 * TDirectoryInformation    directory information type
 */
template<typename TKLocations>
int fm_map(IOConfig const & io_cfg, PrimerConfig & primer_cfg, TKLocations & locations)
{
    std::string const s1 = io_cfg.get_index_dir().string();
    std::string const s2 = io_cfg.get_mapping_dir().string();
    locations.clear();
    TLocations loc_per_K;
    using TKLocationsKey = typename TKLocations::key_type;
    using TKLocationsValue = typename TKLocations::mapped_type;
    for (auto K = PRIMER_MIN_LEN; K <= PRIMER_MAX_LEN; ++K)
    {
        std::cout << "STATUS: run genmap::mappability with E = " << primer_cfg.get_error() << std::endl;
        std::cout << "INFO: K = " << K << std::endl;
        // Remark: csv flag triggers `csvComputation` and therefore the population of the (TLocations) locations vector!
        char const * argv[11] = {"map",
            "-I", s1.c_str(), "-O", s2.c_str(),
            "-K", std::to_string(K).c_str(),
            "-E", std::to_string(primer_cfg.get_error()).c_str(),
            "--csvRAM", "-fl"};

        mappabilityMain<TLocations>(11, argv, loc_per_K, primer_cfg.get_digamma());
        typename TKLocations::iterator it_hint = locations.begin();
        // inserting map pair using hint
        for (auto it = loc_per_K.begin(); it != loc_per_K.end(); ++it)
        {
            TKLocationsKey const key = std::make_tuple(seqan::getValueI1<TSeqNo, TSeqPos>(it->first), seqan::getValueI2<TSeqNo, TSeqPos>(it->first), K);
            TKLocationsValue const value = it->second;
            it_hint = locations.insert(it_hint, std::pair<TKLocationsKey, TKLocationsValue>(key, value));
        }
    }
    return 0;
}

} // namespace priset
