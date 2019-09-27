#include <array>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <experimental/filesystem>
#include <fstream>
#include <numeric>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

#include "priset.hpp"
#include "types.hpp"

namespace fs = std::experimental::filesystem;

// g++ ../PriSeT/src/complete_test.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -lstdc++fs -Wno-unknown-pragmas -o complete_test

struct setup
{
    std::string lib_dir = fs::canonical("../PriSeT/src/tests/library/131221").string();
    //std::cout << "lib_dir = " << lib_dir << std::endl;
    std::string work_dir = fs::canonical("../PriSeT/src/tests/work/131221").string();

    fs::path idx_dir = work_dir + "/index";
    fs::path idx_zip = work_dir + "/index.zip";
    fs::path tmp_dir = work_dir + "/tmp";

    up()
    {
        // unzip index.zip into same named directory
        //std::system(("unzip -n -d " + work_dir + " " + idx_zip.string()).c_str());
        // create tmp dir
        if (!fs::exists(temp_dir) && !fs::create_directory(tmp_dir))
            std::cout << "ERROR: could not create tmp_dir = " << tmp_dir << std::endl;
        std::cout << "lib_dir in setup = " << lib_dir << std::endl;
        std::cout << "work_dir in setup = " << work_dir << std::endl;

    }

    void down()
    {
        // delete index dir
        if (fs::remove_all(idx_dir))
            std::cout << "ERROR: could not remove idx_dir = " << idx_dir << std::endl;
        // delete tmp dir
        if (fs::remove_all(tmp_dir))
            std::cout << "ERROR: could not remove tmp_dir = " << tmp_dir << std::endl;
    }
};

/* Measure runtime for PriSeT components */
void run()
{
    setup su{};

    unsigned const argc = 5;
    char * const argv[argc] = {"priset", "-l", &su.lib_dir[0], "-w", &su.work_dir[0]};
    for (unsigned i = 0; i < argc; ++i) std::cout << argv[i] << " ";
    std::cout << std::endl;

    priset_main(argc, argv);
    std::cout << "MESSAGE: ... done." << std::endl;

}


int main(/*int argc, char ** argv*/)
{
    run();

    return 0;
}
