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

int run_genmap(io_config & cfg)
{
    // check if source directory exists and contains one *.fa and one *.tax file
    DIR *dirp;
    struct dirent *dp;
    std::cout << cfg.fasta_file << std::endl;
    if ((dirp = opendir(cfg.fasta_file.c_str())) == NULL)
        std::cout << "ERROR: " << SRC_DIR_ERROR << std::endl, exit(0);
    bool fasta_set = false, tax_set = false;
    do
    {
        errno = 0;
        if ((dp = readdir(dirp)) != NULL)
        {
            if (!fasta_set && dp->d_reclen > 3 && strstr(dp->d_name, cfg.suffix_fasta))
                cfg.fasta_file += "/" + std::string(dp->d_name), fasta_set = true;
            else if (!tax_set && dp->d_reclen > 4 && strstr(dp->d_name, cfg.suffix_tax))
                cfg.tax_file += "/" + std::string(dp->d_name), tax_set = true;
        }
        else
        {
            if (errno == 0) {
                closedir(dirp);
                std::cout << "ERROR: " << NO_SRC_ERROR << std::endl, exit(0);
            }
            closedir(dirp);
            std::cout << "ERROR: " << SRC_READ_ERROR << std::endl, exit(0);
        }
    }
    while (dp != NULL && !fasta_set && !tax_set);

    std::cout << "Found fasta file:\t" << cfg.fasta_file << "\nFound taxonomy file:\t" << cfg.tax_file << std::endl;

    // create working directory, note that index and map subdirectories have to non-existent
    cfg.genmap_idx_dir += "/genmap_idx";
    cfg.genmap_map_dir += "/genmap_map";

    std::cout << "Create working directory: \t" << cfg.work_dir << std::endl;
    // system(const char*)
    char cmd_rm[50], cmd_mkdir[50];
    sprintf(cmd_rm, "rm -r %s", cfg.work_dir.c_str());
    sprintf(cmd_mkdir, "mkdir -p %s", cfg.work_dir.c_str());
    std::cout << cmd_rm << std::endl;
    std::cout << cmd_mkdir << std::endl;

    if (system(cmd_rm) || system(cmd_mkdir))
        std::cout << "ERROR: " << WRK_DIR_ERROR << std::endl, exit(0);

    std::cout << "Set genmap index directory: " << cfg.genmap_idx_dir << "\nSet genmap mapping directory: " << cfg.genmap_map_dir << std::endl;

    // create FM index and store in cfg.genmap_idx_dir
    // suppress string to char * conversion warning with -Wno-write-strings compiler flag
    pid_t pid;
    if ((pid = fork()) == -1)
        std::cout << "ERROR: " << FORK_ERROR << std::endl, exit(0);
    if (pid == 0) {
        execl(cfg.genmap_bin.c_str(), "genmap", "index", "-F", &cfg.fasta_file[0u], "-I", &cfg.genmap_idx_dir[0u], NULL);
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
            system(&("mkdir " + cfg.genmap_map_dir)[0u]);
            // create mapping of k-mers
            execl(cfg.genmap_bin.c_str(), "genmap", "map", "-I", &cfg.genmap_idx_dir[0u], "-O", &(cfg.genmap_map_dir)[0u], "-K", "18", "-E", "1", "-r", "-d", "-fl", NULL);
            std::cout << "ERROR: " << EXECV_ERROR << std::endl, exit(0);
        }
        else
            wait(NULL);
    }
    return 0;
}
