#include <iostream>
#include <experimental/filesystem>
#include <fstream>
#include <regex>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

#include "filter.hpp"
#include "gui.hpp"
#include "primer_cfg_type.hpp"
#include "types.hpp"
#include "utilities.hpp"

namespace fs = std::experimental::filesystem;

// g++ ../PriSeT/src/test.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -lstdc++fs -o test

struct setup
{
    // TODO: make this runnable with arbitrarily located build folders
    fs::path const lib_dir = "../PriSet/src/tests/library";
    fs::path const work_dir = "../PriSet/src/tests/work";
    priset::io_cfg_type io_cfg;
    priset::primer_cfg_type primer_cfg;
    priset::TKmerLocations kmer_locations;
    priset::TKmerMap kmer_map;

    setup() : io_cfg{lib_dir, work_dir}, primer_cfg{}
    {
        // init primer configurator
        // k1: [(1,2), (1,75)], i.e. kmer1 occurs in  sequence 1 at position 2 and 75
        // k2: [(1,5), (1,80)], i.e. kmer2 occurs in sequence 2 at positions 20 and 80
        // -k1[2]-k2[20]-------------k1[75]/k2[80]
        primer_cfg.set_primer_length_range(priset::primer_cfg_type::size_interval_type{4, 8});
        primer_cfg.set_transcript_range(priset::primer_cfg_type::size_interval_type{50,800});

        priset::TLocation loc1_kmer1{1, 2};
        priset::TLocation loc2_kmer1{1, 75};
        priset::TLocation loc1_kmer2{1, 5};
        priset::TLocation loc2_kmer2{1, 80};

        // init locations
        priset::TKmerLocation::TLocationVec loc_vec1{{loc1_kmer1, loc2_kmer1}};
        priset::TKmerLocation::TLocationVec loc_vec2{{loc1_kmer2, loc2_kmer2}};
        kmer_locations.push_back(priset::TKmerLocation(1, loc_vec1));
        kmer_locations.push_back(priset::TKmerLocation(2, loc_vec2));

        kmer_map[1] = priset::TKmer{1, "AAAA", 12.0};
        kmer_map[2] = priset::TKmer{2, "ACCC", 13.0};
    }

};

void gui_test()
{
    setup up{};
    if (priset::gui::generate_app(up.io_cfg) && priset::gui::compile_app(up.io_cfg))
        std::cout << "success" << std::endl;
    else
        std::cout << "failed to generate and compile app" << std::endl;
}

void combine_test()
{
    setup s{};

    priset::TKmerPairs pairs{};
    priset::combine(s.primer_cfg, s.kmer_locations, s.kmer_map, pairs);
    priset::print_pairs(pairs, s.kmer_map);
}

void create_table_test()
{
    fs::path const lib_dir = "../PriSet/src//tests/library";
    std::cout << "lib_dir = " << lib_dir << " exists: " << fs::exists(lib_dir) << std::endl;
    std::cout << fs::canonical(lib_dir) << std::endl;
    setup s{};

    priset::TKmerPairs pairs{};
    create_table(s.io_cfg, s.kmer_locations, s.kmer_map, pairs);
}

int main()
{
    gui_test();
}
