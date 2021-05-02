// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================
//          Author: Marie Hoffmann <ozymandiaz147 AT gmail.com>
//          Manual: https://github.com/mariehoffmann/PriSeT

#pragma once

#include <experimental/filesystem>
#include <iostream>
#include <unistd.h>

#include "types/IOConfig.hpp"
#include "types/PrimerConfig.hpp"

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
    std::string usage_string = "Usage: %s [-i|s] -l <dir_library> -w <dir_work> [-K <word_length>] [-E <errors>]\n";

    using size_type = PrimerConfig::size_type;

    // Parsed library directory.
    std::string lib_dir;

    // Parsed working directory.
    std::string work_dir;

    // Compute index only, do not compute mappability.
    bool idx_only{0};

    // Parsed flag indicating if index re-computation shall be omitted.
    bool skip_idx{0};

    // Number of errors allowed for k-mer.
    PrimerConfig::size_type E{0};

    // Flags for initializing io configurator.
    bool flag_lib{0}, flag_work{0};

    // Flags for initializing primer configurator.
    bool flag_E{0};
    //
    void parse_arguments(unsigned argc, char * const * argv, PrimerConfig & primer_cfg, IOConfig & io_cfg)
    {
        int opt;
        for (unsigned i = 0; i < argc; ++i) std::cout << argv[i] << " ";
        std::cout << std::endl;

        // TODO: continue here, use only i,s,e as flags
        // l (lib_dir), w (work_dir), i (index only), s (skip_idx), E (error), K (kmer length), colon indicates argument
        while ((opt = getopt(argc, argv, "l:w:isE:")) != -1)
        {
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
                    skip_idx = 0;
                    break;
                case 's':
                    skip_idx = 1;
                    break;
                case 'E':
                    flag_E = 1;
                    E = atoi(optarg);
                    break;
                default: /* '?' */
                    std::cout << "unknown argument opt = " << opt << std::endl;
                    fprintf(stderr, &usage_string[0], argv[0]), exit(EXIT_FAILURE);
             }
        }
        if (!(flag_lib && flag_work))
        {
            fprintf(stderr, &usage_string[0], argv[0]), exit(EXIT_FAILURE);
        }
        // init IO configurator
        io_cfg.assign(lib_dir, work_dir, skip_idx);
        flag_E ? primer_cfg.set_error(E) : (void) (NULL);
    }

public:

    options() = default;
    options(int argc, char * const * argv, PrimerConfig & primer_cfg, IOConfig & io_cfg)
    {
        parse_arguments(argc, argv, primer_cfg, io_cfg);
        // library successfully read, set library and clade sizes
        primer_cfg.set_library_size(io_cfg.get_library_size());
        primer_cfg.set_species_count(io_cfg.get_species_count());
    }

};
}  // namespace priset
