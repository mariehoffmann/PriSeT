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

#include "algorithm.hpp"
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
    Solver(IOConfig & _io_cfg, PrimerConfig & _primer_cfg) :
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

    // maximal score
    size_t C_max{0};              

    // Map of sorted result group indices. key = size, value = index in groups list.
    std::map<size_t, std::vector<size_t>> group_idcs_srtd;

    // State variable for looping over results, either groups or sorted groups.
    using State = std::tuple<size_t, decltype(std::crbegin(group_idcs_srtd)), size_t>;

    // State variable for iterating over results, either groups or sorted groups.
    State state{0, std::crbegin(group_idcs_srtd), 0};

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
        {
            auto key = groups[i].get_taxa_count();
            if (group_idcs_srtd.count(key))
                group_idcs_srtd.at(key).push_back(i);
            else
                group_idcs_srtd[key] = std::vector<size_t>{i};
        }
            // group_idcs_srtd[groups[i].get_taxa_count()] = i;
        is_srtd_by_coverage = true;
        is_srtd_by_frequency = false;
    }

    // Sort result groups by their frequency w.r.t. distinct sequences they match.
    // TODO: group_icds_srtd key not unique!
    void sort_results_by_frequency()
    {
        group_idcs_srtd.clear();
        for (size_t i = 0; i < groups.size(); ++i)
        {
            auto key = groups[i].get_sequence_count();
            if (group_idcs_srtd.count(key))
                group_idcs_srtd.at(key).push_back(i);
            else
                group_idcs_srtd[key] = std::vector<size_t>{i};
        }
            // group_idcs_srtd[groups[i].get_sequence_count()] = i;
        is_srtd_by_coverage = false;
        is_srtd_by_frequency = true;
    }

    // Get result output header for csv output
    const std::string get_header() const noexcept
    {
        return "#name, forward (5'-3'), Tm_fwd, CG_fwd, reverse (5'-3'), Tm_rev, \
        CG_rev, coverage, frequency, [taxa]\n";
    }
    
    /*
    * Iterate over groups (iterator stored in state.first) or sorted groups 
    * (iterator stored in state.second).
    */
    std::pair<bool, Group<TSeqNoMap>> get_next_result()
    {
        auto [ idx1, idx2, pos ] = state;
        // TODO: write unit test, iterator seems not be forwarded
        // if (is_srtd_by_coverage || is_srtd_by_frequency)
        // {
        //     if (idx2 == group_idcs_srtd.rend())
        //         return std::pair<bool, Group<TSeqNoMap>>{false, Group<TSeqNoMap>()};
        //     size_t group_idx = (idx2->second).at(pos);
        //     std::pair<bool, Group<TSeqNoMap>> p{true, groups.at(group_idx)};
        //     if (pos == (idx2->second).size() - 1)
        //         state = std::make_tuple(idx1, ++idx2, 0);
        //     else
        //         state = std::make_tuple(idx1, idx2, ++pos);
        //     return p;
        // } // else
        if (idx1 == groups.size())
            return std::pair<bool, Group<TSeqNoMap>>{false, Group<TSeqNoMap>()};
        std::pair<bool, Group<TSeqNoMap>> p{true, groups.at(idx1)};
        state = std::make_tuple(++idx1, idx2, 0);
        return p;
    }


    std::string generate_statistics()
    {
        std::string info = "Frequency Step,Transform and Filter Step,Combine Step,Pair Filter Step\n";
        info += std::to_string(kmer_counts[0]) + "," + std::to_string(kmer_counts[1]) + ",";
        info += std::to_string(kmer_counts[2]) + "," + std::to_string(kmer_counts[3]) + ",";
        return info;
    }

    // as_string: return also std::string if true else omit
    std::string generate_table(bool as_string=false)
    {
        std::string s{""};
        state = State{0, std::crbegin(group_idcs_srtd), 0};
        auto [success, group] = get_next_result();
        std::ofstream ofs;
        ofs.open(io_cfg.get_result_file());
        // TODO: use get_next_result to iterate over optionally sorted groups
        // while (success)
        // {
        //     s += group.to_string();
        //     ofs << group.to_string();
        //     auto [success, group] = get_next_result();
        // }
        int head = 20;
        if (!groups.size())
        {
            std::cout << "INFO: no primer pairs found with the current settings" << std::endl;
            return s;
        }
        ofs << groups.at(0).get_header();
        std::string aux;
        for (auto group : groups)
        {
            aux = group.to_string();
            if (head > 0)
            {
                s += aux;
                --head;
            }
            ofs << aux;
        }
        ofs.close();
        std::cout << "INFO: output written to " << io_cfg.get_result_file() << std::endl;
        return s;
    }

    void run()
    {
        // 1. Index computation
        if (!io_cfg.get_skip_idx())
        {
            fm_index(io_cfg);
            return;
        }

        // 2. k-mer frequency computation
        fm_map<TKLocations>(io_cfg, primer_cfg, locations);

        // 3. Transform and filter pairs.
        TReferences references;
        TKmerIDs kmerIDs;

        transform_and_filter<TKLocations, TSeqNoMap, TKmerIDs>(io_cfg, primer_cfg, locations, references, seqNo_map, kmerIDs, taxa_by_seqNo_cx, kmer_counts);
        uint64_t tkmer_ctr{0};
        for (auto kmerID_list : kmerIDs)
            tkmer_ctr += kmerID_list.size();

        std::cout << "STATUS\tfound " << tkmer_ctr << " TKmerID\n";
        
        // 4. Combine frequent k-mers to form pairs reference-wise.
        // Container type for storing Pairs.
        using PrimerPairList = std::vector<PrimerPair>;
        PrimerPairList pairs;
        combine<PrimerPairList, TKmerIDs>(references, kmerIDs, pairs, kmer_counts);
        std::cout << "STATUS\tfound " << pairs.size() << " TKmerID pairs\n";
        // 5. Filter and unpack
        filter_and_unpack_pairs<TSeqNoMap, TKmerIDs, PrimerPairList, PrimerPairUnpackedList>(
            io_cfg, 
            primer_cfg, 
            seqNo_map, 
            kmerIDs, 
            pairs, 
            pairs_unpacked
        );
        std::cout << "STATUS\t " << pairs_unpacked.size() << " unpacked and filtered pairs remain" << std::endl; 

        // 6. Write primer pairs to memory
        std::cout << generate_table(true);
    }

    bool as_groups()
    {
        if (!pairs_unpacked.size())
            return false;

        groups.clear();
        for (auto p : pairs_unpacked)
        {
            Group<TSeqNoMap> g(&io_cfg, &seqNo_map, p);
            groups.push_back(g);
        }
        return true;
    }

    // greedy grouping by clustering first primer_set_size primers sorted by coverage.
    virtual bool group_by_max_coverage_greedy()
    {
        if (!pairs_unpacked.size())
            return false;
        groups.clear();
        sort(pairs_unpacked.begin(), pairs_unpacked.end(),
        [](PrimerPairUnpacked<TSeqNoMap> & p, PrimerPairUnpacked<TSeqNoMap> & q)
            {return p.get_species_count() < q.get_species_count();});
        for (auto it = std::crbegin(pairs_unpacked); it != std::crend(pairs_unpacked) - primer_cfg.get_primer_set_size() + 1; ++it)
        {
            PrimerPairUnpackedList ppu{it, it + primer_cfg.get_primer_set_size()};
            Group<TSeqNoMap> group{&io_cfg, &seqNo_map, ppu};
            groups.push_back(group);
        }
        return true;
    }

    // greedy grouping by clustering first primer_set_size primers sorted by coverage.
    bool group_by_max_frequency()
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

    // Solve exactly for the 100 highest covering primer pairs.
    // Rank by largest combined coverage.
    bool group_by_max_coverage_exact()
    {
        // determine rank of species
        std::unordered_map<Taxid, size_t> taxid2rank;
        Taxid taxid;
        size_t rank{0}; // position in common species set
        while ((taxid = io_cfg.get_next_species()))
            taxid2rank[taxid] = rank++;
        std::map<size_t, uint8_t> hitcount2idx;
        // cache best 100 primer pairs in terms of independent coverage
        std::vector<std::vector<bool>> sp_top;
        for (PrimerPairUnpacked<TSeqNoMap> ppu : pairs_unpacked)
        {
            Taxid taxid;
            std::vector<bool> hits(taxid2rank.size(), 0);
            size_t hits_ctr{0};
            while ((taxid = ppu.get_next_species()))
            {
                hits[taxid2rank.at(taxid)] = 1;
                ++hits_ctr;
            }
            if (hitcount2idx.size() < 100) // add
            {
                sp_top.push_back(hits);
                hitcount2idx[hits_ctr] = sp_top.size() - 1;
            }
            else // possibly replace vector with lowest hit count
            {
                std::map<size_t, uint8_t>::iterator it = hitcount2idx.begin();
                if (it->first < hits_ctr)
                {
                    uint8_t idx = it->second;
                    sp_top[idx] = hits;
                    hitcount2idx.erase(it);
                    hitcount2idx[hits_ctr] = idx;
                }
            }
        }
        return (groups.size()) ? true : false;
    }

    void sort_groups_by_coverage()
    {
        assert(false); // TODO: enable upon provisioning of taxid file
        std::sort(groups.begin(), groups.end(), [](Group<TSeqNoMap> & g, Group<TSeqNoMap> & h)
            {return g.get_taxa_count() < h.get_taxa_count();});
    }

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
