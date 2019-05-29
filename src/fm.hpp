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
template<typename io_cfg_type>
int fm_index(io_cfg_type & io_cfg)
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
 * io_cfg_type              I/O configurator type
 * primer_cfg_type          primer configurator type
 * TLocations               type for storing locations
 * TDirectoryInformation    directory information type
 * fasta_header_type        container type for storing fasta header lines
 * fasta_length_type        container type for storing fasta entry lengths (for txt.concat)
 */
template<typename TsequenceNames, typename TsequenceLengths>
int fm_map(io_cfg_type & io_cfg, primer_cfg_type & primer_cfg, TLocations & locations, TDirectoryInformation & directoryInformation, TsequenceNames & sequenceNames, TsequenceLengths & sequenceLengths)
{
    // omit file I/O
    using key_type = typename TLocations::key_type;
    using TSeqNo = typename seqan::Value<key_type, 1>::Type;
    using size_interval_type = typename primer_cfg_type::size_interval_type;
    // seqan::Alloc - for direct memory mapping use seqan::MMap<> instad of seqan::Alloc<>, relevant for benchmarking,
    // since loading index of human genome from disk to main memory may take several minutes
    TFMIndexConfig::SAMPLING = 10;
    // load index
    TIndex index;

    // load index
    if (!genmap::detail::open(index, seqan::toCString(std::string(io_cfg.get_index_base_path())), seqan::OPEN_RDONLY))
        std::cout << "Error in loading index to index obj.\n", exit(0);

    // set directory information
    if (!seqan::open(directoryInformation, seqan::toCString(std::string(io_cfg.get_index_base_path_ids())), seqan::OPEN_RDONLY))
        std::cout << "Error in loading index.ids to directoryInformation obj.\n", exit(0);
    seqan::appendValue(directoryInformation, "dummy.entry;0;chromosomename"); // dummy entry enforces that the mappability is
    uint16_t ctr = 0; // continue here: what is stored in dirInfo
    for (auto it = begin(directoryInformation); it != end(directoryInformation); ++it, ++ctr)
        std::cout << ctr << ": " << (*it) << std::endl;

    // remains empty when excludePseudo == false (store fileIDs for counting matches once per file!)
    std::vector<TSeqNo> mappingSeqIdFile(length(directoryInformation) - 1);
    std::string fastaFile = std::get<0>(retrieveDirectoryInformationLine(directoryInformation[0]));
    std::cout << "fastaFile = " << fastaFile << std::endl;

    // set search parameters (see genmap/src/common.hpp), set in mappabilityMain
    uint16_t K = primer_cfg.template get_primer_length_range<size_interval_type>().first;
    std::cout << "K set to " << K;
    // TODO: when allowing errors, overlap size is different, see formular in mappability.hpp
    uint16_t overlap = std::ceil(0.3 * float(K));
    uint8_t threads = 1; // TODO: set to omp_get_max_threads();
    bool revCompl = false;
    bool excludePseudo = false;
    SearchParams searchParams = SearchParams{K, overlap, threads, revCompl, excludePseudo};

    bool mmap = false;
    bool indels = false;
    bool wigFile = false; // group files into mergable flags, i.e., BED | WIG, etc.
    bool bedFile = false;
    bool rawFile = true;
    bool txtFile = true;
    bool csvFile = true;
    OutputType outputType = OutputType::mappability;
    bool directory = false;
    bool verbose = false;
    CharString indexPath = seqan::CharString(std::string(io_cfg.get_index_base_path()));
    CharString outputPath = ".";
    CharString alphabet = "dna5";
    uint32_t seqNoWidth = 16;
    uint32_t maxSeqLengthWidth = 32;
    uint32_t totalLengthWidth = 32;
    unsigned errors = 0;
    unsigned sampling = 10;

    Options opt{mmap, indels, wigFile, bedFile, rawFile, txtFile, csvFile, outputType, directory, verbose, indexPath, outputPath, alphabet, seqNoWidth, maxSeqLengthWidth, totalLengthWidth, errors, sampling};
    // flag for indicating that index is built on entire directory
    //bool directoryFlag = false;

    // open index file directory
    std::cout << "io_cfg.index_dir : " << io_cfg.get_index_dir() << std::endl;

    run1<TLocations, seqan::Dna5>(locations, opt, searchParams);


    /*
    for (uint64_t i = 0; i < seqan::length(directoryInformation); ++i)  // line 218
    {
        auto const row = retrieveDirectoryInformationLine(directoryInformation[i]);
        std::cout << "row = " << std::get<0>(row) << ", " << std::get<1>(row) << std::endl;
        if (std::get<0>(row) != fastaFile)
        {
            std::cout << std::get<0>(row) << " != fastaFile\n";
            // fastaInfix is now 'text' in mappability.hpp!
            auto const & fastaInfix = seqan::infixWithLength(text.concat, startPos, fastaFileLength);
            std::cout << "fastaInfix = " << fastaInfix << ", fastaFileLength = " << fastaFileLength << std::endl;
            double start = get_wall_time();
            // run<TDistance, value_type, csvComputation, TSeqNo, TSeqPos>(index, fastaInfix, opt, searchParams, fastaFile, sequenceNames, sequenceLengths, directoryInformation, mappingSeqIdFile);
            // TODO: continue here with text -> fastaInfix (see renaming in caller cascade in mappability.hpp)
            std::vector<TValue> c(length(fastaInfix), 0); // line 143, filled by computeMappability, purpose?

            // set <errors=0, csvComputation=false>()
            //run<seqan::HammingDistance, uint16_t, true, TSeqNo, TSeqPos>(index, fastaInfix, opt, searchParams, fastaFile, chromosomeNames, chromosomeLengths, directoryInformation, mappingSeqIdFile);

            //computeMappability<0, false>(index, fastaInfix, c, searchParams, directoryFlag, sequenceLengths, locations, mappingSeqIdFile);
            std::cout << "Mappability computed in " << (round((get_wall_time() - start) * 100.0) / 100.0) << " seconds\n";

            startPos += fastaFileLength;
            fastaFile = std::get<0>(row);
            fastaFileLength = 0;
            for (uint32_t i = 0; i < length(sequenceNames); ++i)
            {
                unsigned id = positionToId(sequenceNames, i);
                std::cout << "chromosomeName = " << valueById(sequenceNames, id);
                id = positionToId(sequenceLengths, i);
                std::cout << ", chromosomeLength = " << valueById(sequenceLengths, id);

            }

            clear(sequenceNames);
            clear(sequenceLengths);
        }
        std::cout << "append len = " << std::get<1>(row) << " to fastaFileLength\n";
        // accumulate fasta entry lengths
        fastaFileLength += std::get<1>(row);
        std::cout << "fastaFileLength acc = " << fastaFileLength << std::endl;
        std::cout << "append chromosomeName = " << std::get<2>(row) << ", chromLength = " << std::get<1>(row) << std::endl;
        appendValue(sequenceNames, std::get<2>(row));
        appendValue(sequenceLengths, std::get<1>(row));
    }*/

    return 0;
}

} // namespace priset
