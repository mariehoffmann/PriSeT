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

    fs::path idx_dir = work_dir + "index";
    fs::path idx_zip = work_dir + "index.zip";
    fs::path tmp_dir = work_dir + "tmp";

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
    uint8_t const K1 = 16;
    uint8_t const K2 = 16;
    std::array<std::array<size_t, priset::TIMEIT::SIZE>, K2 - K1 + 1> runtimes;
    uint8_t const argc = 8;
    for (auto k = K1; k <= K2; ++k)
    {
        std::cout << "workdir as string: " << su.work_dir << std::endl;
        char * const argv[argc] = {"priset", "-l", &su.lib_dir[0], "-w", &su.work_dir[0], "-K", &std::to_string(k)[0], "-s"};

        std::cout << "MESSAGE: start run with K = " << k << std::endl;

        priset_main(argc, argv, &runtimes[k - K1]);
        std::cout << "MESSAGE: ... done." << k << std::endl;
    }

    std::cout << "K\tMAP\tTRANSFORM\tFILTER1\tCOMBINER\tFILTER2\t|\tSUM" << std::string(40, '-') << "\n";
    for (auto k = K1; k <= K2; ++k)
        std::cout << k << '\t' << runtimes[k-K1][priset::TIMEIT::MAP] << '\t' <<
            '\t' << runtimes[k-K1][priset::TIMEIT::TRANSFORM] << '\t' << runtimes[k-K1][priset::TIMEIT::FILTER1] <<
            '\t' << runtimes[k-K1][priset::TIMEIT::COMBINER] << '\t' << runtimes[k-K1][priset::TIMEIT::FILTER2] <<
            "\t|\t" << std::accumulate(std::cbegin(runtimes[k-K1]), std::cend(runtimes[k-K1]), 0) << '\n';
}


int main(/*int argc, char ** argv*/)
{
    timeit();

    return 0;
}
