#pragma once

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <experimental/filesystem>
#include <fstream>
#include <numeric>
#include <regex>
#include <sstream>
#include <string>
#include <sys/wait.h>
#include <tuple>
#include <unistd.h>
#include <unordered_map>
#include <vector>

// #include "argument_parser.hpp"
#include "filter.hpp"
#include "fm.hpp"
#include "types/all.hpp"
#include "solver_fast.hpp"
#include "solver_brute_force.hpp"
#include "utilities.hpp"

namespace fs = std::experimental::filesystem;

using namespace priset;

struct Solver
{
private:

    // Stores for each reference the encoded kmers in order of occurrence.
    using TKmerIDs = std::vector<std::deque<TKmerID>>;

    // A k-mer location augmented by the information about K.
    using TKLocation = std::tuple<priset::TSeqNo, priset::TSeqPos, priset::TKmerLength>;

    // The map of k-mer locations augmented by information of k value to preserve
    // key uniqueness. Contains mappings for all values of K in range (see
    // PrimerConfig.hpp).
    using TKLocations = std::map<TKLocation,
            std::pair<std::vector<TLocation>, std::vector<TLocation>>>;

    // Container type for storing results.
    using ResultList = std::vector<Result>;

    // Translates sequences identifiers (seqNo) in use to a contiguous range (seqNo_cx).
    // Dictionary is bidirectional: seqNo -> seqNo_cx and inverse add a leading one
    // to the compressed key: (1 << 63 | seqNo_cx) -> seqNo.
    // Background: some sequences produce no k-mers and therefore no space should be
    // reserved in its bit transformation.
    using TSeqNoMap = std::unordered_map<TSeqNo, TSeqNo>;

    // Reference to io configurator.
    IOConfig & io_cfg;

    // Reference to primer configurator.
    PrimerConfig & primer_cfg;

    vector<size_t> & primerIDs;

    // vector<vector<bool>> & solutions;  // optimal primer combination
    size_t C_max{0};              // maximal score

    PairList pairs;               // list of type Pair

    // Result collector.
    std::vector<std::vector<Result>> & solutions;

    // Map of sorted solution indices. key = size, value = solution index.
    std::map<size_t, size_t> solutions_srtd;

    // State variable for looping over results, either solutions or sorted solutions.
    std::pair<size_t, std::map<size_t, size_t>> state{0, solutions.rbegin()};

    // K-Mer location map.
    TKLocations locations;

    // TDirectoryInformation directoryInformation;

    // The FASTA header lines of the reference sequences.
    TSequenceNames sequenceNames;

    // The lengths of the reference sequences.
    TSequenceLengths sequenceLengths;

public:
    Solver(IOConfig & _io_cfg, PrimerConfig & _primer_cfg) :
    io_cfg(_io_cfg), primer_cfg(_primer_cfg) {}

    virtual solve();

    void set_solutions(std::vector<ResultList> const & _solutions)
    {
        solutions = _solutions;
    }

    void sort_results_by_coverage()
    {
        for (size_t id = 0; id < solutions.size(); ++id)
        {
            std::unordered_set<uint64_t> taxon_set;
            for (Result result : solutions[id])
            {
                std::vector<uint64_t> taxa = result.get_taxa();
                taxon_set.insert(taxa.begin(), taxa.end());
            }
            solutions_srtd[taxon_set.size()] = id;
        }
    }

    //
    void sort_results_by_frequency()
    {
        for (size_t id = 0; id < solutions.size(); ++id)
        {
            std::unordered_set<uint64_t> ref_set;
            for (Result result : solutions[id])
            {
                std::vector<uint64_t> refs = result.get_references();
                taxa_set.insert(refs.begin(), refs.end());
            }
            solutions_srtd[ref_set.size()] = id;
        }
    }

    // Get result output header for csv output
    constexpr std::string get_header() const noexcept
    {
        return "#name, forward (5'-3'), Tm_fwd, CG_fwd, reverse (5'-3'), Tm_rev, \
        CG_rev, coverage, frequency, [taxa]\n";
    }

    std::pair<bool, ResultList> get_next_result()
    {
        if (solutions_srtd.size())
        {
            if (state.second == solutions_srtd.rend())
                return pair<bool, ResultList>{false, ResultList(0)};
            return pair<bool, ResultList>{true, solutions.at(state.second++)};
        }
        if (state == results.size())
            return pair<bool, ResultList>{false, ResultList(0)};
        return pair<bool, ResultList>{true, results.at(state.first++)};
        }
    }

    bool generate_table()
    {
        state{0, solutions.rbegin()};
        auto [success, result] = get_next_result();
        ofstream ofs;
        ofs.open(io_cfg.get_result_file());
        while (success)
        {
            ofs << result.to_string();
            auto [success, result] = get_next_result();
        }
        ofs.close();
        std::cout << "INFO: output written to " << io_cfg.get_result_file() << std::endl;
        return true;
    }

private:
    // copy app code here
    run()
    {
        // 1. Optional index computation
        if (!io_cfg.skip_idx_flag)
            fm_index(io_cfg);

        // 2. k-mer frequency computation
        fm_map(io_cfg, primer_cfg, locations);

        // 3. Transform and filter pairs.
        TReferences references;
        TKmerIDs kmerIDs;
        TSeqNoMap seqNoMap;
        transform_and_filter(io_cfg, locations, references, seqNoMap, kmerIDs);

        // 4. Combine frequent pairs reference-wise.
        using PairList = PairList<Pair<CombinePattern<TKmerID, TKmerLength>>>;
        PairList pairs;
        combine<PairList>(references, kmerIDs, pairs);

        // 5. Maximize for coverage
        using PairLists = std::vector<PairList>;
        PairLists pairs_grouped;
        optimize_coverage<PairList, PairLists>(pairs, primer_cfg.get_primer_s, pairs_grouped);

        // 6. Filter pairs.
        ResultList results;

        // TODO: rewrite  filter_pairs to work with pair_freqs
        filter_and_retransform<PairLists, ResultList>(references, kmerIDs, pairs_grouped, results);
    }

};
