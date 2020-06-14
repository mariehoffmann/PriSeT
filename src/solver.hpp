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

#include "filter.hpp"
#include "fm.hpp"
#include "types/all.hpp"
#include "utilities.hpp"

namespace fs = std::experimental::filesystem;

namespace priset
{

struct Solver
{
    // Reference to io configurator.
    IOConfig & io_cfg;

    // Reference to primer configurator.
    PrimerConfig & primer_cfg;

    // Constructor with member variable setting.
    Solver(IOConfig & _io_cfg, PrimerConfig _primer_cfg) :
        io_cfg(_io_cfg), primer_cfg(_primer_cfg) {}

    // Stores for each reference the encoded kmers in order of occurrence.
    using TKmerIDs = std::vector<std::deque<TKmerID>>;

    // A k-mer location augmented by the information about K.
    using TKLocation = std::tuple<priset::TSeqNo, priset::TSeqPos, priset::TKmerLength>;

    // The map of k-mer locations augmented by information of k value to preserve
    // key uniqueness. Contains mappings for all values of K in range (see
    // PrimerConfig.hpp).
    using TKLocations = std::map<TKLocation,
            std::pair<std::vector<TLocation>, std::vector<TLocation>>>;

    // Translates sequences identifiers (seqNo) in use to a contiguous range (seqNo_cx).
    using TSeqNoMap = std::unordered_map<TSeqNo, TSeqNo>;

    // The container type for groups. template<typename TSeqNoMap>
    using Groups = std::vector<Group<TSeqNoMap>>;

    // The container type for unpacked pairs.
    using PrimerPairUnpackedList = std::vector<PrimerPairUnpacked<TSeqNoMap>>;


private:
    // K-Mer location map.
    TKLocations locations;

    // Dictionary is bidirectional: seqNo -> seqNo_cx and inverse add a leading one
    // to the compressed key: (1 << 63 | seqNo_cx) -> seqNo.
    // Background: some sequences produce no k-mers and therefore no space should be
    // reserved in its bit transformation.
    TSeqNoMap seqNo_map;

    // The container for groups.
    Groups groups;

    // The container for PrimerPairUnpacked.
    PrimerPairUnpackedList pairs_unpacked;

    // Stores for each sequence its assigned taxon. The index corresponds to the
    // continuous range of sequence identifiers accessable via the seqNo_map.
    std::vector<Taxid> taxa_by_seqNo_cx;

    std::vector<size_t> & primerIDs();

    size_t C_max{0};              // maximal score

    // Map of sorted result group indices. key = size, value = index in groups list.
    std::map<size_t, size_t> group_idcs_srtd;

    // State variable for looping over results, either groups or sorted groups.
    using State = std::pair<size_t, decltype(std::crbegin(group_idcs_srtd))>;

    // State variable for iterating over results.
    State state{0, std::crbegin(group_idcs_srtd)};

    // Flag indicating if indices in group_idcs_srtd represent sorting by coverage.
    bool is_srtd_by_coverage{false};

    // Flag indicating if indices in group_idcs_srtd represent sorting by frequency.
    bool is_srtd_by_frequency{false};

    // K-mer counts at each step for analysis purposes.
    uint64_t kmer_counts[4] = {};

    // TDirectoryInformation directoryInformation;

    // The FASTA header lines of the reference sequences.
    TSequenceNames sequenceNames;

    // The lengths of the reference sequences.
    TSequenceLengths sequenceLengths;

public:

    virtual void solve()
    {
        std::cerr << "Missing implementation of void solve() in your derived class!" << std::endl;
    };

    template<typename TSeqNoMap>
    void add_group(Group<TSeqNoMap> const & group) noexcept
    {
        groups.push_back(group);
    }

    // Sort result groups by taxonomic coverage.
    void sort_results_by_coverage()
    {
        group_idcs_srtd.clear();
        for (size_t i = 0; i < groups.size(); ++i)
            group_idcs_srtd[groups[i].get_taxa_count()] = i;
        is_srtd_by_coverage = true;
        is_srtd_by_frequency = false;
    }

    // Sort result groups by their frequency w.r.t. distinct sequences they match.
    void sort_results_by_frequency()
    {
        group_idcs_srtd.clear();
        for (size_t i = 0; i < groups.size(); ++i)
            group_idcs_srtd[groups[i].get_sequence_count()] = i;
        is_srtd_by_coverage = false;
        is_srtd_by_frequency = true;
    }

    // Get result output header for csv output
    const std::string get_header() const noexcept
    {
        return "#name, forward (5'-3'), Tm_fwd, CG_fwd, reverse (5'-3'), Tm_rev, \
        CG_rev, coverage, frequency, [taxa]\n";
    }

    std::pair<bool, Group<TSeqNoMap>> get_next_result()
    {
        if (group_idcs_srtd.size())
        {
            if (state.second == group_idcs_srtd.rend())
                return std::pair<bool, Group<TSeqNoMap>>{false, Group<TSeqNoMap>()};
            std::pair<bool, Group<TSeqNoMap>> p{true, groups.at((state.second++)->second)};
            // ++state.second;
            return p;
        }
        if (state.first == groups.size())
            return std::pair<bool, Group<TSeqNoMap>>{false, Group<TSeqNoMap>()};
        return std::pair<bool, Group<TSeqNoMap>>{true, groups.at(state.first++)};
    }

    std::string generate_statistics()
    {
        std::string info = "Frequency Step,Transform and Filter Step,Combine Step,Pair Filter Step\n";
        info += std::to_string(kmer_counts[0]) + "," + std::to_string(kmer_counts[1]) + ",";
        info += std::to_string(kmer_counts[2]) + "," + std::to_string(kmer_counts[3]) + ",";
        return info;
    }

    bool generate_table()
    {
        state = State{0, std::crbegin(group_idcs_srtd)};
        auto [success, group] = get_next_result();
        std::ofstream ofs;
        ofs.open(io_cfg.get_result_file());
        while (success)
        {
            ofs << group.to_string();
            auto [success, group] = get_next_result();
        }
        ofs.close();
        std::cout << "INFO: output written to " << io_cfg.get_result_file() << std::endl;
        return true;
    }

    // copy app code here
    void run()
    {
        // 1. Optional index computation
        if (!io_cfg.skip_FM_idx())
            fm_index(io_cfg);

        // 2. k-mer frequency computation
        fm_map<TKLocations>(io_cfg, primer_cfg, locations);

        // 3. Transform and filter pairs.
        TReferences references;
        TKmerIDs kmerIDs;

        transform_and_filter<TKLocations, TSeqNoMap, TKmerIDs>(io_cfg, primer_cfg, locations, references, seqNo_map, kmerIDs, taxa_by_seqNo_cx, kmer_counts);

        // 4. Combine frequent k-mers to form pairs reference-wise.
        // Container type for storing Pairs.
        using PrimerPairList = std::vector<PrimerPair>;
        PrimerPairList pairs;
        combine<PrimerPairList, TKmerIDs>(references, kmerIDs, pairs, kmer_counts);

        // 5. Filter and unpack
        filter_and_unpack_pairs<TSeqNoMap, TKmerIDs, PrimerPairList, PrimerPairUnpackedList>(io_cfg, primer_cfg, seqNo_map, kmerIDs, pairs, pairs_unpacked);

    }

    bool as_groups()
    {
        if (!pairs_unpacked.size())
            return false;
        std::transform(pairs_unpacked.begin(), pairs_unpacked.end(), groups.begin(),
            [&](PrimerPairUnpacked<TSeqNoMap> & p)
            {return Group<TSeqNoMap>(&io_cfg, &seqNo_map, p);});
        return true;
    }

    // greedy grouping by clustering first primer_set_size primers sorted by coverage.
    virtual bool group_by_max_coverage()
    {
        if (!pairs_unpacked.size())
            return false;
        groups.clear();
        sort(pairs_unpacked.begin(), pairs_unpacked.end(),
        [](PrimerPairUnpacked<TSeqNoMap> & p, PrimerPairUnpacked<TSeqNoMap> & q)
            {return p.get_coverage() < q.get_coverage();});
        for (auto it = std::crbegin(pairs_unpacked); it != std::crend(pairs_unpacked) - primer_cfg.get_primer_set_size() + 1; ++it)
        {
            PrimerPairUnpackedList ppu{it, it + primer_cfg.get_primer_set_size()};
            Group<TSeqNoMap> group{&io_cfg, &seqNo_map, ppu};
            groups.push_back(group);
        }

        return true;
    }

    // greedy grouping by clustering first primer_set_size primers sorted by coverage.
    virtual bool group_by_max_frequency()
    {
        if (!pairs_unpacked.size())
            return false;
        groups.clear();
        sort(pairs_unpacked.begin(), pairs_unpacked.end(),
        [](PrimerPairUnpacked<TSeqNoMap> & p, PrimerPairUnpacked<TSeqNoMap> & q)
            {return p.get_frequency() < q.get_frequency();});
        for (auto it = std::crbegin(pairs_unpacked); it != std::crend(pairs_unpacked) - primer_cfg.get_primer_set_size() + 1; ++it)
        {
            PrimerPairUnpackedList ppu{it, it + primer_cfg.get_primer_set_size()};
            groups.push_back(Group<TSeqNoMap>{&io_cfg, &seqNo_map, ppu});
        }

        return true;
    }

    // // TODO: continue here
    // // 6. Maximize for coverage
    // using Groups = std::vector<Group>;
    // Groups groups;
    // optimize_coverage<PrimerPairList, Groups>(io_cfg, primer_cfg, pairs_unpacked, groups, kmer_counts);
    //
    // // 7. Filter pairs.
    // // TODO: rewrite  filter_pairs to work with pair_freqs
    // using PrimerPairGroupsUnpack =
    // filter_groups<TKmerIDs, PrimerPairGroups, PrimerPairGroupsUnpack>(primer_cfg, references, kmerIDs, pairs_grouped, pairs_unpacked_grouped, kmer_counts);
    //
    // // 8. convert groups containing KMerID encoded pairs into DNA sequences
    // Group results;
    // retransform<Group>(pairs_unpacked_grouped, results);

};

} // namespace priset
