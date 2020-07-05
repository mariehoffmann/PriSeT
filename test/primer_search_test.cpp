#include <array>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <experimental/filesystem>
#include <fstream>
#include <numeric>
#include <regex>
#include <sstream>
#include <string>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

#include "../src/argument_parser.hpp"
#include "../src/filter.hpp"
#include "../src/fm.hpp"
#include "../src/IOConfig.hpp"
#include "../src/PrimerConfig.hpp"
#include "../src/priset.hpp"
#include "../src/types.hpp"
#include "../src/utilities.hpp"

//g++ ../PriSeT/tests/primer_search_test.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -lstdc++fs -Wno-unknown-pragmas -lstdc++fs -DNDEBUG -O3 -I/Users/troja/include -L/Users/troja/lib -lsdsl -ldivsufsort -o primer_search_test

namespace fs = std::experimental::filesystem;

template<typename TMatches>
void init_primer_search(fs::path & primer_file, TMatches & matches)
{
    std::ifstream infile(primer_file);
    std::string line;
    std::pair<std::string, std::string> key;
    while (std::getline(infile, line))
    {
        std::cout << "line = " << line << std::endl;
        if (line.compare(0, 1, ">"))
        {
            std::cout << "ERROR: expect line starting with '>'\n";
            break;
        }
        auto header = line;
        if (!std::getline(infile, line))
        {
            std::cout << "ERROR: unexpected end of file\n";
            break;
        }
        std::cout << "next line = " << line << std::endl;
        key = std::pair<std::string, std::string>{header.substr(1, header.size() - 1), line};
        std::cout << "key = " << key.first << ", " << key.second << std::endl;
        matches[key] = 0;
    }
}

int main()
{
    std::unordered_set<std::string> kmers_unique;
    kmers_unique.insert("ACTTTCGTTCTTGATTGA");   // 1
    kmers_unique.insert("AAACTTGAAGAGGCG");      // 0
    kmers_unique.insert("AAACTTGAAGTGGCGG");     // 1
    kmers_unique.insert("AAAGTCTCGTGAGTGC");     // 1
    kmers_unique.insert("TAAAGTCTCGTGAGTGC");    // 0
    kmers_unique.insert("GGACAGAAAGACCCTATG");   // 1
    kmers_unique.insert("AGATCAGGCTGTTATCC");    // 0
    kmers_unique.insert("GGTCAACAAATCATAAAGATATTGG"); // 1
    kmers_unique.insert("GCTTGTCTCAAAGATTAAGCC"); // 1
    kmers_unique.insert("GCCTGCTGCCTTCCTTGGA");  // 1
    kmers_unique.insert("AAACTCAAAGAGACGG");     // 1
    kmers_unique.insert("AAACTCGAAGAGGCGG");     // 1
    kmers_unique.insert("AAACTTAAAGAGACGG");     // 1
    kmers_unique.insert("AAACTTAAAGTGACGG");     // 1
    kmers_unique.insert("AAACTTGAAGTGGCGG");     // 1

    fs::path primers_known = "/Users/troja/git/PriSet_git2/PriSeT/tests/ecology_primers_regex.fasta";
    if (!fs::exists(primers_known))
    {
        std::cout << "ERROR: file path incorrect '" << primers_known << "'\n";
        exit(-1);
    }
    std::map<std::pair<std::string, std::string>, uint16_t> matches;
    init_primer_search(primers_known, matches);
    std::cout << "matches.size() = " << matches.size() << std::endl;
    // Iterate over regex and unique kmer strings and report matches
    for (auto it = matches.begin(); it != matches.end(); ++it)
    {
        for (std::string const & kmer_str : kmers_unique)
        {
            std::regex const primer_rx(it->first.second);
            bool found = std::regex_match(kmer_str, primer_rx);
            if (found)
                ++(it->second);
        }
        std::cout << "(" << it->first.first << ", " << it->first.second << "): " << it->second << std::endl;
    }
    return 0;
}
