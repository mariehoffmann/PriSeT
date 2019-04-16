#include <iostream>
#include <fstream>
#include <regex>
#include <vector>

// g++ test.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -o test

int main()
{
    std::string tax_file = "/Users/troja/tmp/priset/src/taxonomy.tax";
    std::ifstream ifs;
    ifs.open(tax_file, std::ifstream::in);
    std::cout << "ifs is open: " << ifs.is_open() << std::endl;
    // store flat, i.e., [pid1, cid1, pid2, cid2, ...]
    std::vector<std::string> raw_taxa;
    std::string cell = "";
    // Iterate through each line and split the content using delimeter
    while (getline(ifs, cell, ','))
    {
        std::cout << cell << ", ";
        raw_taxa.push_back(cell);
    }
    std::cout << std::endl;
    std::cout << "raw_taxa.size() = " << raw_taxa.size() << std::endl;
    ifs.close();

    std::regex suffix_fasta_rx("\\.(fa|fasta)");
}
