#pragma once

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <iostream>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

#include "io_config.hpp"


/*
# Build genmap binary used for FM index building
cmake ../genmap -DCMAKE_BUILD_TYPE=Release -DGENMAP_NATIVE_BUILD=OFF
make genmap
*/

/*
 * Create FM index with `genmap` and store in io_cfg.genmap_idx_dir
 */
int fm_index(io_config & io_cfg, primer_config & primer_cfg)
{

    
    pid_t pid;
    if ((pid = fork()) == -1)
        std::cout << "ERROR: " << FORK_ERROR << std::endl, exit(0);
    if (pid == 0) {
        execl(io_cfg.genmap_bin.c_str(), "genmap", "index", "-F", &io_cfg.fasta_file[0u], "-I", &io_cfg.genmap_idx_dir[0u], NULL);
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
 */
int fm_map(io_config & io_cfg, primer_config & primer_cfg)
{
    // TODO: allow optional file I/O
    //io_cfg.genmap_map_dir += "/genmap_map";
    pid_t pid;

    /*
    // ./bin/genmap map -I ~/tmp/genmap -O ~/tmp/genmap -K 8 -E 1 --raw -t -d -fl
    if ((pid = fork()) == -1)
        std::cout << "ERROR: " << FORK_ERROR << std::endl, exit(0);
    if (pid == 0) {
        std::cout << "start mapping ...\n";
        // create mapping output directory, which has to be already existent for genmap!
        system(&("mkdir " + io_cfg.genmap_map_dir)[0u]);
        // create mapping of k-mers

        execl(io_cfg.genmap_bin.c_str(), "genmap", "map",
            "-I", &io_cfg.genmap_idx_dir[0u],
            "-O", &(cfg.genmap_map_dir)[0u],
            "-K", std::string(primer_cfg.get_primer_length_range.min),
            "-E", "1", "-r", "-d", "-fl", NULL);
        std::cout << "ERROR: " << EXECV_ERROR << std::endl, exit(0);
    }
    else
        wait(NULL);
    */
    return 0;
}
