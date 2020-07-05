#include <array>
#include <chrono>
// #include <cstdlib>
#include <iostream>
#include <experimental/filesystem>
// #include <fstream>
// #include <numeric>
#include <sys/wait.h>
#include <unistd.h>
#include <unordered_map>
#include <vector>

#include <seqan/basic.h>

#include "../src/argument_parser.hpp"
#include "../src/errors.hpp"
#include "../src/filter.hpp"
#include "../src/fm.hpp"
#include "../src/gui.hpp"
#include "../src/IOConfig.hpp"
#include "../src/output.hpp"
#include "../src/PrimerConfig.hpp"
#include "../src/taxonomy.hpp"
#include "../src/types.hpp"
#include "../src/utilities.hpp"

namespace fs = std::experimental::filesystem;
using namespace priset;

// static constexpr bool outputProgress = false;

// -mpopcnt gives speedup of 5, otherwise call __popcountdi2 called for __builtin_popcountll
// g++ ../PriSeT/tests/performance_test.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -lstdc++fs -Wno-unknown-pragmas -lstdc++fs -DNDEBUG -O3 -mpopcnt -I/Users/troja/include -L/Users/troja/lib -lsdsl -ldivsufsort -o performance_test


/* Measure runtime for PriSeT components */
int main(/*int argc, char ** argv*/)
{
    priset::TKmerID kmerID = 6917529143072389695;
    std::cout << "kmerID = " << kmerID2str(kmerID) << std::endl;
    return 0;
}
