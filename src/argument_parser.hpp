// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================
//          Author: Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
//          Manual: https://github.com/mariehoffmann/PriSeT
#pragma once

#include <experimental/filesystem>
#include <iostream>
#include <unistd.h>

#include "io_cfg_type.hpp"
#include "primer_cfg_type.hpp"

namespace fs = std::experimental::filesystem;

namespace priset
{
/*
 * Type to store optional arguments given via command line. Will launch setup of primer and io configurators.
 * Library and working directory info goes to io, and K, E to primer configurator.
 */
struct options
{
private:
    std::string usage_string = "Usage: %s -l <dir_library> -w <dir_work> [-K <word_length>] [-E <errors>]\n";

    using size_type = primer_cfg_type::size_type;

    // Parsed library directory.
    std::string lib_dir;

    // Parsed working directory.
    std::string work_dir;

    // Compute index only, do not compute mappability.
    bool idx_only{0};

    // Parsed flag indicating if index re-computation shall be omitted.
    bool skip_idx{0};

    // Length of k-mer. Per default the lower bound from the primer configurator is taken.
    primer_cfg_type::size_type K;

    // Number of errors allowed for k-mer.
    primer_cfg_type::size_type E{0};

    // Flags for initializing io configurator.
    bool flag_lib{0}, flag_work{0};
    // Flags for initializing primer configurator.
    bool flag_E{0}, flag_K{0};
    //
    void parse_arguments(int argc, char * argv[], primer_cfg_type & primer_cfg, io_cfg_type & io_cfg)
    {
        int opt;

        // l (lib_dir), w (work_dir), i (index only), s (skip_idx), E (error), K (kmer length), colon indicates argument
        while ((opt = getopt(argc, argv, "l:w:isE:K:")) != -1)
        {
            std::cout << "opt = " << opt << std::endl;
            switch (opt)
            {
                case 'l':
                    flag_lib = 1;
                    lib_dir.assign(std::string(optarg));
                    break;
                case 'w':
                    flag_work = 1;
                    work_dir.assign(std::string(optarg));
                    break;
                case 'i':
                    idx_only = 1;
                    break;
                case 's':
                    skip_idx = 1;
                    break;
                case 'E':
                    flag_E = 1;
                    E = atoi(optarg);
                    break;
                case 'K':
                    flag_K = 1;
                    K = atoi(optarg);
                    break;
                default: /* '?' */
                    std::cout << "unknown argument opt = " << opt << std::endl;
                    fprintf(stderr, &usage_string[0], argv[0]), exit(EXIT_FAILURE);
             }
        }
        if (!(flag_lib && flag_work))
            fprintf(stderr, &usage_string[0], argv[0]), exit(EXIT_FAILURE);
        // init io configurator
        io_cfg.assign(lib_dir, work_dir, idx_only, skip_idx);
        flag_K ? primer_cfg.set_primer_length_range(K) : (void) (NULL);
        flag_E ? primer_cfg.set_error(E) : (void) (NULL);
    }

public:

    options() = default;
    options(int argc, char * argv[], primer_cfg_type & primer_cfg, io_cfg_type & io_cfg)
    {
        parse_arguments(argc, argv, primer_cfg, io_cfg);
    }

};
}  // namespace priset
