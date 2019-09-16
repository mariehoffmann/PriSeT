#include <sdsl/bit_vectors.hpp>

#include <iostream>

using namespace std;
using namespace sdsl;

//g++ ../PriSeT/src/test2.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -lstdc++fs -DNDEBUG -O3 -I/Users/troja/include -L/Users/troja/lib -lsdsl -ldivsufsort -o test2

int main()
{

    for (uint64_t k = 1; k <= 0; ++k)

    {
        std::cout << "loop" << std::endl;
    }
    std::vector<std::pair<uint64_t, uint64_t>> pairs;
    bit_vector b(10000, 0);
    b[0] = 1;
    b[5] = 1;
    b[60] = 1;
    b[65] = 1;
    b[300] = 1;

    sdsl::rank_support_v5<1, 1> r1s(&b);
    sdsl::select_support_mcl<1> s1s(&b);

    cout << "i\tb(i)\trank(i)\tselect(rank(i))\n";
    for (uint64_t i = 0; i < 10; ++i)
        cout <<  i << "\t" << b[i] << "\t" << r1s.rank(i)  <<  "\t"  << s1s.select(1) << endl;

    cout << "select(1) = " << s1s.select(1) << std::endl;
    cout << "select(2) = " << s1s.select(2) << std::endl;
    cout << "select(3) = " << s1s.select(3) << std::endl;
    cout << "select(4) = " << s1s.select(4) << std::endl;
    cout << "select(10000) = " 

    uint64_t offset_max = 500;
    uint64_t offset_min = 40;
    uint64_t K = 16;
    for (uint64_t r = 1; r < r1s.rank(b.size()); ++r)
    {
        uint64_t idx = s1s.select(r); // position of r-th k-mer
        // min offset
        cout << "idx = " << idx << endl;
        uint64_t w_begin = std::min(b.size()-1, idx + K + offset_min + 1);
        uint64_t w_end = std::min(b.size(), idx + K + offset_max + 1);
        cout << "w_begin = " << w_begin << ", w_end = " << w_end << ", rank(w_begin) = " << r1s.rank(w_begin) << ", rank(w_end) = " << r1s.rank(w_end) << endl;
        cout << "number of ones in window: " << (r1s.rank(w_end) - r1s.rank(w_begin)) << std::endl;
        cout << "combine kmer idx = " << r;
        for (uint64_t k = 1; k <= r1s.rank(w_end) - r1s.rank(w_begin); ++k)
        {
            cout << "with kmer idx2 = " << r1s.rank(w_begin) + k << endl;
            uint64_t idx2 = s1s.select(r1s.rank(w_begin) + k);
            pairs.push_back(std::make_pair(idx, idx2));
            std::cout << "pair: (" << idx << ", " << idx2 << ")\n";
        }
    }
    return 0;
}
