#include <cassert>
#include <chrono>
#include <iostream>
#include <experimental/filesystem>
#include <fstream>
#include <regex>
#include <sys/wait.h>
#include <unistd.h>
#include <unordered_map>
#include <vector>

#include <seqan/basic.h>
#include <seqan/file.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

#include "../src/filter.hpp"
#include "../src/PrimerConfig.hpp"
#include "../src/types.hpp"
#include "../src/utilities.hpp"

using namespace priset;

// g++ ../PriSeT/tests/types_test.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -lstdc++fs -DNDEBUG -O3 -I/Users/troja/include -L/Users/troja/lib -lsdsl -ldivsufsort -I .. -o types_test

void test_CombinePattern()
{
    CombinePattern cp{};

    std::vector<std::pair<uint8_t, uint8_t>> cs;
    for (uint8_t i = 0; i < PREFIX_SIZE; ++i)
    {
        for (uint8_t j = 0; j < PREFIX_SIZE; ++j)
        {
            cp.set(ONE_LSHIFT_63 >> i, ONE_LSHIFT_63 >> j);
            cp.get_combinations(cs);
            if (!cs.size())
                std::cout << "ERROR: combination for (" << int(i) <<  ", " << int(j) << "not set\n";
            else if (cs.size() > 1)
                std::cout << "ERROR: cs too large with " << cs.size() << ", expect 0." << std::endl;
            else
            {
                if (cs.size() == 1 && cs.begin()->first == i && cs.begin()->second == j)
                {
                    std::cout << "SUCCESS for set\n";

                    // unset with
                    cp.reset(i, j);
                    cp.get_combinations(cs);
                    if (!cs.size())
                        std::cout << "SUCCESS for reset\n";
                    else
                        std::cout << "ERROR: expect reset combinations\n";
                }
                else
                    std::cout << "ERROR: expect (" << int(i) << ", " << int(j) << "), got (" << int(cs.begin()->first) << ", " << int(cs.begin()->second) << ")" << std::endl;
            }
        }
    }
}

int main()
{
    test_CombinePattern();
    return 0;
}
