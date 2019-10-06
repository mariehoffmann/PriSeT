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

#include "../src/priset.hpp"
#include "../src/types.hpp"

namespace fs = std::experimental::filesystem;

// g++ ../PriSeT/tests/performance_test.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -lstdc++fs -Wno-unknown-pragmas -lstdc++fs -DNDEBUG -O3 -I/Users/troja/include -L/Users/troja/lib -lsdsl -ldivsufsort -o performance_test

struct setup
{
    std::string lib_dir = fs::canonical("../PriSeT/tests/library/3041").string();
    //std::cout << "lib_dir = " << lib_dir << std::endl;
    std::string work_dir = fs::canonical("../PriSeT/tests/work/3041").string();

    fs::path idx_dir = work_dir + "/index";
    fs::path idx_zip = work_dir + "/index.zip";
    fs::path tmp_dir = work_dir + "/tmp";

    setup()
    {
        // unzip index.zip into same named directory
        std::system(("unzip -n -d " + work_dir + " " + idx_zip.string()).c_str());
        // clear and create tmp dir
        std::cout << "tmpdir exists: " << fs::exists(tmp_dir) << std::endl;
        if (fs::exists(tmp_dir))
        {
            std::cout << "tmp_dir exists, delete it\n";
            fs::remove_all(tmp_dir);
        }
        if (!fs::create_directory(tmp_dir))
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

    std::cout << "K\tMAP\t\tFILTER1_TRANSFORM\tCOMBINE_FILTER2\tPAIR_FREQ\t|\tSUM [μs]\n" << std::string(100, '_') << "\n";
    std::cout << "[" << 16 << ":" << 25 << "]\t" << runtimes[priset::TIMEIT::MAP] << "\t" <<
            '\t' << runtimes[priset::TIMEIT::FILTER1_TRANSFORM] <<
            '\t' << runtimes[priset::TIMEIT::COMBINE_FILTER2] << '\t' << runtimes[priset::TIMEIT::PAIR_FREQ] <<
            "\t|\t" << std::accumulate(std::cbegin(runtimes), std::cend(runtimes), 0) << '\n';

}


int main(/*int argc, char ** argv*/)
{
    timeit();

    return 0;
}