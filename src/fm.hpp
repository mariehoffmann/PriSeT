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

#include <seqan/arg_parse.h>
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
 * io_config        I/O configurator type
 * primer_config    primer configurator type
 * TSeqNo           genmap's sequence counter type
 * TSeqPos          genmap's sequence position type
 */
template<typename io_config, typename primer_config, typename TLocations>
int fm_map(io_config & io_cfg, primer_config & primer_cfg, TLocations & locations)
{
    // omit file I/O
    // ./bin/genmap map -I ~/tmp/genmap -O ~/tmp/genmap -K 8 -E 1 --raw -t -d -fl
/*
    int argc = 12; // Todo: add or remove program name
    const char* argv[12] = { "map",
        "-I", &std::string(io_cfg.get_index_dir())[0u],         // index directory
        "-O", &std::string(io_cfg.get_mapping_dir())[0u],       // output directory
        "-K", &std::to_string(primer_cfg.get_primer_length_range().min)[0u],  // k-mer length, default 18 bp
        "-E", "0",      // number of errors
        "-r",           // raw output, i.e. binary of std::vector<uint16_t>
        "-d",           // detailed csv (omit later)
        "-fl"           // frequency-large, i.e. uint16_t (max = (1<<16)-1)
    };
    */// Value<TTuple, I>::Type;
    using key_type = typename TLocations::key_type;
    using TSeqNo = typename seqan::Value<key_type, 1>::Type;
    using TSeqPos = typename seqan::Value<key_type, 2>::Type;
    // FM index mapper settings (see mappability.hpp for further information)
    // parameter definitions for setting the index type
    typedef String<seqan::TChar, seqan::TAllocConfig> TString;
    typedef seqan::StringSet<TString, seqan::Owner<seqan::ConcatDirect<SizeSpec_<TSeqNo, TSeqPos> > > > TStringSet;
    using TBWTLen = uint64_t;
    using TFMIndexConfig = TGemMapFastFMIndexConfig<TBWTLen>;
    // set sampling rate, TODO: verify if 10 is ok
    TFMIndexConfig::SAMPLING = 10;
    // set index type
    using TIndex = seqan::Index<TStringSet, TBiIndexConfig<TFMIndexConfig> >;
    // load index
    TIndex index;
    open(index, seqan::toCString(io_cfg.index_dir), OPEN_RDONLY);

    auto const & text = indexText(index);
    std::vector<value_type> c(length(text), 0);
    // set search parameters (SearchParams struct defined in common.hpp)
    unsigned length = primer_cfg.get_primer_length_range().min;
    unsigned overlap = length * 0.7;
    unsigned threads = 1;
    bool revCompl = false;
    bool excludePseudo = true;
    SearchParams searchParams = SearchParams{length, overlap, threads, revCompl, excludePseudo};

    // flag for indicating that index is built on entire directory
    bool directoryFlag = false;
    // only needed for file I/O (except for csv), TODO: verify
    seqan::StringSet<uint64_t> chromosomeLengths;

    // fill directory information struct
    seqan::StringSet<seqan::CharString, seqan::Owner<seqan::ConcatDirect<> > > directoryInformation;
    // open index file  ".ids"
    open(directoryInformation, seqan::toCString(std::string(io_cfg.get_index_dir() + ".ids"), seqan::OPEN_RDONLY);
    appendValue(directoryInformation, "dummy.entry;0;chromosomename");
    std::vector<TSeqNo> mappingSeqIdFile(length(directoryInformation) - 1);
    if (searchParams.excludePseudo)
    {
        uint64_t fastaId = 0;
        std::string fastaFile = std::get<0>(retrieveDirectoryInformationLine(directoryInformation[0]));
        for (uint64_t i = 0; i < length(directoryInformation) - 1; ++i)
        {
            auto const row = retrieveDirectoryInformationLine(directoryInformation[i]);
            if (std::get<0>(row) != fastaFile)
            {
                fastaFile = std::get<0>(row);
                ++fastaId;
            }
            mappingSeqIdFile[i] = fastaId;
        }
    }

    double start = get_wall_time();
    // set <errors=0, csvComputation=true>()
    computeMappability<0, true>(index, text, c, searchParams, directoryFlag, chromLengths, locations, mappingSeqIdFile);
    std::cout << "Mappability computed in " << (round((get_wall_time() - start) * 100.0) / 100.0) << " seconds\n";

    return 0;
}

} // namespace priset
