#include <iostream>
#include <fstream>
#include <regex>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

#include "filter.hpp"
#include "primer_cfg_type.hpp"
#include "types.hpp"
#include "utilities.hpp"

// g++ test.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -o test

// k1: [(1,2), (1,75)], i.e. kmer1 occurs in  sequence 1 at position 2 and 75
// k2: [(1,20), (1,80)], i.e. kmer2 occurs in sequence 2 at positions 20 and 80
// -k1[2]-k2[20]-------------k1[75]/k2[80]
void test_combine()
{
    // init primer configurator
    priset::primer_cfg_type primer_cfg{};
    primer_cfg.set_primer_length_range(priset::primer_cfg_type::size_interval_type{18, 24});
    primer_cfg.set_transcript_range(priset::primer_cfg_type::size_interval_type{50,800});
    // init locations
    priset::TKmerLocations kmer_locations;
    priset::TLocation loc1_kmer1{1, 2};
    priset::TLocation loc2_kmer1{1, 75};
    priset::TLocation loc1_kmer2{1, 20};
    priset::TLocation loc2_kmer2{1, 80};

    using TLocationVec = typename std::vector<priset::TLocation>;
    kmer_locations.push_back(std::make_pair<priset::TKmerID, TLocationVec>(1, TLocationVec{loc1_kmer1, loc2_kmer1}));
    kmer_locations.push_back(std::make_pair<priset::TKmerID, TLocationVec>(2, TLocationVec{loc1_kmer2, loc2_kmer2}));

    priset::TKmerMap kmer_map{{1, priset::TKmer{1, "AAAA", 12.0}}, {2, priset::TKmer{2, "CCCC", 13.0}}};
    priset::TKmerPairs pairs{};
    priset::combine(primer_cfg, kmer_locations, kmer_map, pairs);
    priset::print_pairs(pairs, kmer_map);
}


int main()
{
    test_combine();

}
