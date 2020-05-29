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
#include "../src/io_cfg_type.hpp"
#include "../src/primer_cfg_type.hpp"
#include "../src/priset.hpp"
#include "../src/types.hpp"
#include "../src/utilities.hpp"

namespace fs = std::experimental::filesystem;

// g++ ../PriSeT/tests/clade_X_test.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -lstdc++fs -Wno-unknown-pragmas -lstdc++fs -DNDEBUG -O3 -I/Users/troja/include -L/Users/troja/lib -lsdsl -ldivsufsort


void load_primers_known(fs::path const & primer_file, std::vector<TPrimerPair> & primer_pairs)
{
    primer_pairs.push_back(TPrimerPair{"DIV4 fwd", dna_encoder("GCGGTAATTCCAGCTCCAATAG"), "DIV4 rev", dna_encoder("TATTCGTATTCCATTGTCAGAG")});
    primer_pairs.push_back(TPrimerPair{"EUK14 fwd", dna_encoder("CAGCAGCCGCGGTAATTCC"), "EUK14 rev", dna_encoder("GCTTAATTTGACTCAACACGGG")});
    primer_pairs.push_back(TPrimerPair{"nSSU fwd", dna_encoder("GCTTGTCTCAAAGATTAAGCC"), "nSSU rev", dna_encoder("TCCAAGGAAGGCAGCAGGC")});
    primer_pairs.push_back(TPrimerPair{"V9 fwd", dna_encoder("GTACACACCGCCCGTC"), "V9 rev", dna_encoder("GTAGGTGAACCTGCAGAAGGATCA")});
    primer_pairs.push_back(TPrimerPair{"23S fwd", dna_encoder("GGACAAAAAGACCCTATG"), "23S rev", dna_encoder("GGATAACAGGCTGATCT")});
    primer_pairs.push_back(TPrimerPair{"23S fwd", dna_encoder("GGACAGAAAGACCCTATG"), "23S rev", dna_encoder("GGATAACAGGCTGATCT")});
    // combine unfolded ones
    std::vector<std::string> EUK15_fwds{"CCAGCACCCGCGGTAATTCC", "CCAGCACCTGCGGTAATTCC", "CCAGCAGCCGCGGTAATTCC", "CCAGCAGCTGCGGTAATTCC"};
    std::vector<std::string> EUK15_revs{"ACTTTCGTTCTTGATCAA", "ACTTTCGTTCTTGATCGA", "ACTTTCGTTCTTGATTAA", "ACTTTCGTTCTTGATTGA"};
    uint16_t i = 1;
    for (std::string EUK15_fwd : EUK15_fwds)
    {
        uint16_t j = 1;
        for (std::string EUK15_rev : EUK15_revs)
            primer_pairs.push_back(TPrimerPair{"EUK15 fwd" + std::to_string(i++), dna_encoder(EUK15_fwd), "EUK15 rev " + std::to_string(j++), dna_encoder(EUK15_rev)});
    }

    std::vector<std::string> UNIV_fwds{"AAACTCAAAGAGACGG", "AAACTCAAAGAGGCGG", "AAACTCAAAGTGACGG", "AAACTCAAAGTGGCGG", "AAACTCGAAGAGACGG", "AAACTCGAAGAGGCGG", "AAACTCGAAGTGACGG", "AAACTCGAAGTGGCGG", "AAACTTAAAGAGACGG", "AAACTTAAAGAGGCGG", "AAACTTAAAGTGACGG", "AAACTTAAAGTGGCGG", "AAACTTGAAGAGACGG", "AAACTTGAAGAGGCGG", "AAACTTGAAGTGACGG", "AAACTTGAAGTGGCGG"};
    std::vector<std::string> UNIV_revs{"GCACACACGAGACTTT", "GCACTCACGAGACTTT", "GTACACACGAGACTTT", "GTACTCACGAGACTTT"};
    i = 1;
    for (std::string UNIV_fwd : UNIV_fwds)
    {
        uint16_t j = 1;
        for (std::string UNIV_rev : UNIV_revs)
            primer_pairs.push_back(TPrimerPair{"UNIV fwd" + std::to_string(i++), dna_encoder(UNIV_fwd), "UNIV rev " + std::to_string(j++), dna_encoder(UNIV_rev)});
    }
}


int main(int argc, char ** argv)
{
    std::vector<TPrimerPair> primer_pairs;
    fs::path primers_known_file = "/Users/troja/git/PriSet_git2/PriSeT/tests/ecology_primers_regex.fasta";
    load_primers_known(primers_known_file, primer_pairs);
    for (auto pp : primer_pairs)
        std::cout << pp.name_fwd << ": " << pp.code_fwd << "\n" << pp.name_rev << ": " << pp.code_rev << "\n\n";
    return 0;
}
