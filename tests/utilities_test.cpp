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

#include "../src/combine_types.hpp"
#include "../src/filter.hpp"
#include "../src/primer_cfg_type.hpp"
#include "../src/types.hpp"
#include "../src/utilities.hpp"

using namespace priset;

// g++ ../PriSeT/tests/utilities_test.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -lstdc++fs -DNDEBUG -O3 -I/Users/troja/include -L/Users/troja/lib -lsdsl -ldivsufsort -I .. -o utilities_test

void get_code_test()
{
    // kmer of length 17 with length bits for l = 16, 17
    std::string s17 = "AAAACCCCAAAACCCCT";
    TKmerID kmerID = (3ULL << 62) + dna_encoder(s17);
    uint64_t code = get_code(kmerID, 1ULL << 63);
    /*if (s17.compare(0, 16, dna_decoder(code, 1ULL << 63)) != 0)
        std::cout << "ERROR: expect " << s17.substr(0, 16) << ", got " << dna_decoder(code, 1ULL << 63) << std::endl;
    else
        std::cout << "SUCCESS\n";
    code = get_code(kmerID, 1ULL << 62);
    if (s17.compare(0, 17, dna_decoder(code, 1ULL << 62)) != 0)
        std::cout << "ERROR: expect " << s17 << ", got " << dna_decoder(code, 1ULL << 62) << std::endl;
    else
        std::cout << "SUCCESS\n";
*/
    code = get_code(kmerID, 0);
    if (code != dna_encoder(s17))
        std::cout << "ERROR: expect " << dna_encoder(s17) << ", got " << code << std::endl;
    else
        std::cout << "SUCCESS\n";

    std::string s = "CAGCAGCCGCGGTAATTCC";
    code = get_code(dna_encoder(s));
    std::cout << "get_code(s) = " << code << std::endl;

    std::string s2 = "CAGCAGCCGCGGTAATTCCG";
    auto code2 = get_code(dna_encoder(s), 1ULL << (63-3));
    if (code != code2)
        std::cout << "ERROR: expect " << code << ", got " << code2 << std::endl;
    else
        std::cout << "SUCCESS\n";
}

int main()
{
    /*
    kmerID =	16070356149012225618
    get_code =	369306287669
    dna_decode =	CCCTTTACCTGGTTGATCC
    */
    std::string s = "CAGCAGCCGCGGTAATTCCC";
    uint64_t code = dna_encoder(s);
    std::cout << "s = \t\t" << s << std::endl;
    std::cout << "code = \t\t" << code << std::endl;
    std::cout << "decode = \t" << dna_decoder(code, 1ULL << 60) << std::endl;

    //get_code_test();
    return 0;
}
