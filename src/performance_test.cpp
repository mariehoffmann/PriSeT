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

// g++ ../PriSeT/src/performance_test.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -lstdc++fs -o performance_test

struct setup
{
    std::string lib_dir = fs::canonical("../PriSeT/src/tests/library/3041").string();
    //std::cout << "lib_dir = " << lib_dir << std::endl;
    std::string work_dir = fs::canonical("../PriSeT/src/tests/work/3041").string();

    fs::path idx_dir = work_dir + "/index";
    fs::path idx_zip = work_dir + "/index.zip";
    fs::path tmp_dir = work_dir + "/tmp";

    setup()
    {
        // unzip index.zip into same named directory
        std::system(("unzip -n -d " + work_dir + " " + idx_zip.string()).c_str());
        // create tmp dir
        if (fs::create_directory(tmp_dir))
            std::cout << "ERROR: could not create tmp_dir = " << tmp_dir << std::endl;
        std::cout << "lib_dir in setup = " << lib_dir << std::endl;
        std::cout << "work_dir in setup = " << work_dir << std::endl;

    }

    void cleanup()
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
void timeit()
{
    setup su{};
    std::array<size_t, priset::TIMEIT::SIZE> runtimes;

    unsigned const argc = 6;
    char * const argv[argc] = {"priset", "-l", &su.lib_dir[0], "-w", &su.work_dir[0], "-s"};
    for (unsigned i = 0; i < argc; ++i) std::cout << argv[i] << " ";
    std::cout << std::endl;

    priset_main(argc, argv, &runtimes);
    std::cout << "MESSAGE: ... done." << std::endl;

    std::cout << "K\tMAP\t\tTRANSFORM\tFILTER1\tCOMBINER\tFILTER2\t|\tSUM [μs]\n" << std::string(100, '_') << "\n";
    std::cout << "[" << 16 << ":" << 25 << "]\t" << runtimes[priset::TIMEIT::MAP] << "\t" <<
            '\t' << runtimes[priset::TIMEIT::TRANSFORM] << '\t' << runtimes[priset::TIMEIT::FILTER1] <<
            '\t' << runtimes[priset::TIMEIT::COMBINER] << '\t' << runtimes[priset::TIMEIT::FILTER2] <<
            "\t|\t" << std::accumulate(std::cbegin(runtimes), std::cend(runtimes), 0) << '\n';
}


int main(/*int argc, char ** argv*/)
{
    timeit();

    return 0;
}
