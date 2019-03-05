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

/*
 * usage        compile g++ run.cpp -o run
 *              ./run <src_dir> <work_dir>
 * e.g.         ./run ~/tmp/priset/src ~/tmp/priset/work
 *
 * src_dir      path to folder containing fasta (*.fa) and taxonomy file (*.tax)
 * work_dir     path to store indices, mappings, annotations, and other results
 */
int main(int argc, char** argv)
{
    if (argc < 3)
        std::cout << "error: " << ARG_ERROR << std::endl, exit(0);

    std::string fa_file = argv[1], tax_file = argv[2], genmap_idx_dir, genmap_map_dir;
    char suffix_fasta[] = ".fa", suffix_tax[] = ".tax";

    std::cout << "src directory: " << fa_file << std::endl;
    // check if source directory exists and contains one *.fa and one *.tax file
    DIR *dirp;
    struct dirent *dp;
    std::cout << argv[1] << std::endl;
    if ((dirp = opendir(argv[1])) == NULL)
        std::cout << "error: " << SRC_DIR_ERROR << std::endl, exit(0);
    bool fasta_set = false, tax_set = false;
    do
    {
        errno = 0;
        if ((dp = readdir(dirp)) != NULL)
        {
            std::cout << "current file = " << dp->d_name << " with len = " << dp->d_reclen << std::endl;
            if (dp->d_reclen > 3 && strstr(dp->d_name, suffix_fasta)) //std::equal(dp->d_name + dp->d_reclen - 3, dp->d_name + dp->d_reclen, ".fa"))
                fa_file += "/" + std::string(dp->d_name), fasta_set = true, std::cout<< "set fa_file = " << fa_file << std::endl;
            else if (dp->d_reclen > 4 && strstr(dp->d_name, suffix_tax))
                tax_file += "/" + std::string(dp->d_name), tax_set = true, std::cout << "set tax_file = " << tax_file << std::endl;;
        }
        else
        {
            if (errno == 0) {
                closedir(dirp);
                std::cout << "error: " << NO_SRC_ERROR << std::endl, exit(0);;
            }
            closedir(dirp);
            std::cout << "error: " << SRC_READ_ERROR << std::endl, exit(0);;
        }
    }
    while (dp != NULL && !fasta_set && !tax_set);

    std::cout << "Found fasta file: " << fa_file << "\nFound taxonomy file: " << tax_file << std::endl;

    // create working directory
    dirp = opendir(argv[2]);
    if (!dirp)
        return -2;
    genmap_idx_dir = std::string(argv[1]) + "/genmap_idx";
    genmap_map_dir = std::string(argv[1]) + "/genmap_map";

    system(std::string("mkdir -p " + genmap_idx_dir).c_str());
    system(std::string("mkdir -p " + genmap_map_dir).c_str());


/*
    vector<char*> argv;
    argv.push_back("genmap");
    argv.push_back("index");
    argv.push_back("-F");
    argv.push_back(argv[0]);
    argv.push_back("-I"); // FM index location
    -I ~/tmp/genmap

    ./bin/genmap index -F ~/gitlab/phd/papers/2017/mueggelsee2014/preparation/MSA_18_species/Ceratium_hirundinella/Ceratium_hirundinella.fa -I ~/tmp/genmap
    execv("genmap", "index", , NULL);
*/
    return 0;
}
