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
 * TLocations       type for storing locations
 */
template<typename io_config, typename primer_config, typename TLocations>
int fm_map(io_config & io_cfg, primer_config & primer_cfg, TLocations & locations)
{
    // omit file I/O
    // ./bin/genmap map -I ~/tmp/genmap -O ~/tmp/genmap -K 8 -E 1 --raw -t -d -fl
    // FM index mapper settings (see mappability.hpp for further information)
    // parameter definitions for setting the index type

    using key_type = typename TLocations::key_type;
    using TSeqNo = typename seqan::Value<key_type, 1>::Type;
    using TSeqPos = typename seqan::Value<key_type, 2>::Type;
    // use frequency_large, other from ./genmap map: frequency_small (uint8_t)
    using TValue = uint16_t;
    //using size_type = typename primer_config::size_type;
    using size_interval_type = typename primer_config::size_interval_type;
    // Dna5 alphabet = {A, C, G, T, N}, letters other than dna4 are casted to 'N'
    // seqan::Alloc - for direct memory mapping use seqan::MMap<> instad of seqan::Alloc<>, relevant for benchmarking,
    // since loading index of human genome from disk to main memory may take several minutes
    typedef String<seqan::Dna5, seqan::Alloc<>> TString;
    typedef seqan::StringSet<TString, seqan::Owner<seqan::ConcatDirect<SizeSpec_<TSeqNo, TSeqPos> > > > TStringSet;
    using TBWTLen = uint64_t;
    using TFMIndexConfig = TGemMapFastFMIndexConfig<TBWTLen>;
    // set sampling rate, TODO: verify if 10 is ok
    TFMIndexConfig::SAMPLING = 10;
    // set index type
    using TIndex = seqan::Index<TStringSet, TBiIndexConfig<TFMIndexConfig> >;
    // path to index plus basename without suffix, i.e. <path_to_index>/<basename>
    std::string index_path_base = std::string(io_cfg.get_index_dir()) + "/index";
    std::string index_path_base_ids = index_path_base + ".ids";
    // load index
    TIndex index;
    open(index, seqan::toCString(index_path_base), OPEN_RDONLY); //open(index, toCString(opt.indexPath), OPEN_RDONLY);

    seqan::StringSet<seqan::CharString, seqan::Owner<seqan::ConcatDirect<> > > directoryInformation;
    seqan::open(directoryInformation, seqan::toCString(index_path_base_ids), seqan::OPEN_RDONLY);
    // dummy entry enforces that the mappability is computed for the last file in the while loop.
    seqan::appendValue(directoryInformation, "dummy.entry;0;chromosomename");

    // remains empty when excludePseudo == false (store fileIDs for counting matches once per file!)
    std::vector<TSeqNo> mappingSeqIdFile(length(directoryInformation) - 1);

    auto const & text = indexText(index); // line 213
    std::vector<TValue> c(length(text), 0); // line 143, filled by computeMappability, purpose?
    seqan::StringSet<seqan::CharString, seqan::Owner<seqan::ConcatDirect<> > > chromosomeNames;
    seqan::StringSet<uint64_t> chromosomeLengths; // filled by computeMappability
    uint64_t startPos = 0;
    uint64_t fastaFileLength = 0;
    std::string fastaFile = std::get<0>(retrieveDirectoryInformationLine(directoryInformation[0]));

    // set search parameters (see common.hpp) with SearchParams{length, overlap, threads, revCompl, excludePseudo}, set in mappabilityMain
    SearchParams searchParams = SearchParams{primer_cfg.template get_primer_length_range<size_interval_type>().first, 12, 1, false, false};
    // flag for indicating that index is built on entire directory
    bool directoryFlag = false;

    // open index file directory
    std::cout << "io_cfg.index_dir : " << io_cfg.get_index_dir() << std::endl;

    for (uint64_t i = 0; i < seqan::length(directoryInformation); ++i)  // line 218
    {
        auto const row = retrieveDirectoryInformationLine(directoryInformation[i]);
        std::cout << "row = " << std::get<0>(row) << ", " << std::get<1>(row) << std::endl;
        if (std::get<0>(row) != fastaFile)
        {
            // fastaInfix is now 'text' in mappability.hpp!
            auto const & fastaInfix = seqan::infixWithLength(text.concat, startPos, fastaFileLength);
            double start = get_wall_time();
            // set <errors=0, csvComputation=true>()
            // computeMappability<0, csvComputation>(index, text, c, searchParams, opt.directory, chromLengths, locations, mappingSeqIdFile);
            // TODO: continue here
            computeMappability<0, true>(index, fastaInfix, c, searchParams, directoryFlag, chromosomeLengths, locations, mappingSeqIdFile);
            std::cout << "Mappability computed in " << (round((get_wall_time() - start) * 100.0) / 100.0) << " seconds\n";

            startPos += fastaFileLength;
            fastaFile = std::get<0>(row);
            fastaFileLength = 0;
            clear(chromosomeNames);
            clear(chromosomeLengths);
        }
        fastaFileLength += std::get<1>(row);
        appendValue(chromosomeNames, std::get<2>(row));
        appendValue(chromosomeLengths, std::get<1>(row));
    }

    return 0;
}

} // namespace priset
