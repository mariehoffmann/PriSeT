#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include <dirent.h>
#include <unistd.h>
#include <vector>
#include <sys/types.h>
#include <sys/wait.h>

#define ARG_ERROR -1
#define SRC_DIR_ERROR -2
#define NO_SRC_ERROR -3
#define SRC_READ_ERROR -4
#define WRK_DIR_ERROR -5
#define FORK_ERROR -6
#define EXECV_ERROR -7

/*
 * usage        g++ run.cpp -Wno-write-strings -o run
 *              ./run <path_genmap_bin> <src_dir> <work_dir>
 * e.g.         ./run $GENMAP ~/tmp/priset/src ~/tmp/priset/work
 *
 * src_dir      path to folder containing fasta (*.fa) and taxonomy file (*.tax)
 * work_dir     path to store indices, mappings, annotations, and other results
 */
int main(int argc, char** argv)
{
    if (argc < 4)
        std::cout << "ERROR: " << ARG_ERROR << std::endl, exit(0);

    // set path prefixes for fasta and taxonomy files
    std::string fa_file = argv[2], tax_file = argv[2];
    char suffix_fasta[] = ".fa", suffix_tax[] = ".tax";
    std::string genmap_idx_dir, genmap_map_dir;

    // check if source directory exists and contains one *.fa and one *.tax file
    DIR *dirp;
    struct dirent *dp;
    std::cout << argv[2] << std::endl;
    if ((dirp = opendir(argv[2])) == NULL)
        std::cout << "ERROR: " << SRC_DIR_ERROR << std::endl, exit(0);
    bool fasta_set = false, tax_set = false;
    do
    {
        errno = 0;
        if ((dp = readdir(dirp)) != NULL)
        {
            if (!fasta_set && dp->d_reclen > 3 && strstr(dp->d_name, suffix_fasta))
                fa_file += "/" + std::string(dp->d_name), fasta_set = true;
            else if (!tax_set && dp->d_reclen > 4 && strstr(dp->d_name, suffix_tax))
                tax_file += "/" + std::string(dp->d_name), tax_set = true;
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

    std::cout << "Found fasta file:\t" << fa_file << "\nFound taxonomy file:\t" << tax_file << std::endl;

    // create working directory, note that index and map subdirectories have to non-existent
    genmap_idx_dir = std::string(argv[3]) + "/genmap_idx";
    genmap_map_dir = std::string(argv[3]) + "/genmap_map";

    std::cout << "Create working directory: \t" << argv[3] << std::endl;
    // system(const char*)
    char cmd_rm[50], cmd_mkdir[50];
    sprintf(cmd_rm, "rm -r %s", argv[3]);
    sprintf(cmd_mkdir, "mkdir -p %s", argv[3]);
    std::cout << cmd_rm << std::endl;
    std::cout << cmd_mkdir << std::endl;

    if (system(cmd_rm) || system(cmd_mkdir))
        std::cout << "ERROR: " << WRK_DIR_ERROR << std::endl, exit(0);

    std::cout << "Set genmap index directory: " << genmap_idx_dir << "\nSet genmap mapping directory: " << genmap_map_dir << std::endl;

    // create FM index and store in genmap_idx_dir
    // suppress string to char * conversion warning with -Wno-write-strings compiler flag
    pid_t pid;
    if ((pid = fork()) == -1)
        std::cout << "ERROR: " << FORK_ERROR << std::endl, exit(0);
    if (pid == 0) {
        int ret = execl(argv[1], "genmap", "index", "-F", &fa_file[0u], "-I", &genmap_idx_dir[0u], NULL);
        std::cout << "ERROR: " << EXECV_ERROR << std::endl, exit(0);
    }
    else
    {
        wait(NULL);
        // create mapping of k-mers
        // ./bin/genmap map -I ~/tmp/genmap -O ~/tmp/genmap -K 8 -E 1 --raw -t -d -fl
        if ((pid = fork()) == -1)
            std::cout << "ERROR: " << FORK_ERROR << std::endl, exit(0);
        if (pid == 0) {
            std::cout << "start mapping ...\n";
            int ret = execl(argv[1], "genmap", "map", "-I", &genmap_idx_dir[0u], "-O", &genmap_map_dir[0u], "-K", "18", "-E", "1", "--raw", "-t", "-d", "-fl", NULL);
            std::cout << "ERROR: " << EXECV_ERROR << std::endl, exit(0);
        }
        else
            wait(NULL);
    }

    return 0;
}
