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

#include "types/IOConfig.hpp"


# Build genmap binary used for FM index building
cmake ../genmap -DCMAKE_BUILD_TYPE=Release -DGENMAP_NATIVE_BUILD=OFF
make genmap


// create FM index and map
int FM_index(IOConfig & io_cfg, PrimerConfig & primer_cfg)
{

    // TODO: move directory creation into io configurator
    // create working directory, note that index and map subdirectories have to non-existent
    io_cfg.genmap_idx_dir += "/genmap_idx";
    io_cfg.genmap_map_dir += "/genmap_map";

    std::cout << "Create working directory: \t" << io_cfg.work_dir << std::endl;
    // system(const char*)
    char cmd_rm[50], cmd_mkdir[50];
    sprintf(cmd_rm, "rm -r %s", io_cfg.work_dir.c_str());
    sprintf(cmd_mkdir, "mkdir -p %s", io_cfg.work_dir.c_str());
    std::cout << cmd_rm << std::endl;
    std::cout << cmd_mkdir << std::endl;

    if (system(cmd_rm) || system(cmd_mkdir))
        std::cout << "ERROR: " << WRK_DIR_ERROR << std::endl, exit(0);

    std::cout << "Set genmap index directory: " << io_cfg.genmap_idx_dir << "\nSet genmap mapping directory: " << io_cfg.genmap_map_dir << std::endl;

    // create FM index and store in io_cfg.genmap_idx_dir
    // suppress string to char * conversion warning with -Wno-write-strings compiler flag
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
    }
    return 0;
}

int FM_map(IOConfig & io_cfg, PrimerConfig & primer_cfg)
{
    return 0;
}
