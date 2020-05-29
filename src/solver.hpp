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
#include <unistd.h>
#include <vector>

#include "argument_parser.hpp"
#include "filter.hpp"
#include "fm.hpp"
#include "io_cfg_type.hpp"
#include "primer_cfg_type.hpp"
#include "solver_fast.hpp"
#include "solver_brute_force.hpp"
#include "types.hpp"
#include "utilities.hpp"

namespace fs = std::experimental::filesystem;

using namespace priset;
// using namespace std;

struct Solver
{
private:

    io_cfg_type & io_cfg;

    // Reference to primer configurator.
    primer_cfg_type & primer_cfg;

    vector<size_t> & primerIDs;

    // vector<vector<bool>> & solutions;  // optimal primer combination
    size_t C_max{0};              // maximal score

    TPairList pairs;               // list of type TPair

    // Result collector.
    std::vector<std::vector<TResult>> & solutions;

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

    using TResultList = std::vector<TResult>;

public:
    Solver(io_cfg_type & _io_cfg, primer_cfg_type & _primer_cfg) :
    io_cfg(_io_cfg), primer_cfg(_primer_cfg) {}

    virtual solve();

    void set_solutions(std::vector<TResultList> const & _solutions)
    {
        solutions = _solutions;
    }

    void sort_results_by_coverage()
    {
        for (size_t id = 0; id < solutions.size(); ++id)
        {
            std::unordered_set<uint64_t> taxon_set;
            for (TResult result : solutions[id])
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
            for (TResult result : solutions[id])
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

    std::pair<bool, TResultList> get_next_result()
    {
        if (solutions_srtd.size())
        {
            if (state.second == solutions_srtd.rend())
                return pair<bool, TResultList>{false, TResultList(0)};
            return pair<bool, TResultList>{true, solutions.at(state.second++)};
        }
        if (state == results.size())
            return pair<bool, TResultList>{false, TResultList(0)};
        return pair<bool, TResultList>{true, results.at(state.first++)};
        }
    }

private:
    // copy app code here
    run()
    {
        // 1. k-mer frequency computation
        fm_map(io_cfg, primer_cfg, locations);

        // 2. Transform and filter pairs.
        TReferences references;
        TKmerIDs kmerIDs;
        TSeqNoMap seqNoMap;
        transform_and_filter(io_cfg, locations, references, seqNoMap, kmerIDs);

        // 3. Combine frequent pairs reference-wise.
        using TPairList = TPairList<TPair<TCombinePattern<TKmerID, TKmerLength>>>;
        TPairList pairs;
        combine<TPairList>(references, kmerIDs, pairs);

        // 4. Maximize for coverage
        using TPairLists = std::vector<TPairList>;
        TPairLists pairs_grouped;
        optimize_coverage<TPairList, TPairLists>(pairs, primer_cfg.get_primer_s, pairs_grouped);

        // 5. Filter pairs.
        TResultList results;

        // TODO: rewrite  filter_pairs to work with pair_freqs
        filter_and_retransform<TPairLists, TResultList>(references, kmerIDs, pairs_grouped, results);
    }

};
