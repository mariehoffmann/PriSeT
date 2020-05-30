#include <cassert>
#include <chrono>
#include <iostream>
#include <experimental/filesystem>
#include <fstream>
#include <regex>
#include <string>
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

namespace fs = std::experimental::filesystem;
using namespace priset;

// g++ ../PriSeT/src/tests/dna_encoder_test.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -lstdc++fs -DNDEBUG -O3 -I/Users/troja/include -L/Users/troja/lib -lsdsl -ldivsufsort -I .. -o dna_encoder_test

void run()
{
    //                 ACGT    ACGT    ACGT    AATT    AATT    AATT
    // 0000000000000001000110110001101100011011000011110000111100001111
    // 0000000000000001000110110001101100011011000011110000111100001111
    // 0000000000000001011011000110110001101100001111000011110000111100
    // 3210987654321098765432109876543210987654321098765432109876543210
    // 19,21,23, encoded length 24
    std::string seq = "ACGTACGTACGTAATTAATTAATT";
    TKmerID kmerID = dna_encoder(seq);
    if (kmerID != 0b0000000000000001000110110001101100011011000011110000111100001111)
        std::cout << "ERROR: expect " << 311278208749327ULL << "got " << kmerID << std::endl;
    // add length information
    kmerID += (ONE_LSHIFT_63 >> 3) + (ONE_LSHIFT_63 >> 5) + (ONE_LSHIFT_63 >> 7);

    std::string res = dna_decoder(kmerID, ONE_LSHIFT_63 >> 3);
    if (!res.compare(0, 19, seq))
    {
        std::cout << "ERROR: expect '" << seq.substr(0, 19) << ", got " << res << std::endl;
        exit(0);
    }
    res = dna_decoder(kmerID, ONE_LSHIFT_63 >> 5);
    if (!res.compare(0, 21, seq))
    {
        std::cout << "ERROR: expect '" << seq.substr(0, 19) << ", got " << res << std::endl;
        exit(0);
    }
    res = dna_decoder(kmerID, ONE_LSHIFT_63 >> 7);
    if (!res.compare(0, 23, seq))
    {
        std::cout << "ERROR: expect '" << seq.substr(0, 19) << ", got " << res << std::endl;
        exit(0);
    }
    std::cout << "SUCCESS\n";
}

int main()
{
    run();
    return 0;
}
