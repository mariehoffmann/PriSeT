// #include <algorithm>
// #include <array>
// #include <chrono>
// #include <cstdlib>
// #include <ctime>
// #include <iostream>
// #include <experimental/filesystem>
// #include <fstream>
// #include <numeric>
// #include <regex>
// #include <sstream>
// #include <string>
// #include <sys/wait.h>
// #include <unistd.h>
// #include <vector>
//
// #include "../src/argument_parser.hpp"
// #include "../src/algorithm.hpp"
// #include "../src/fm.hpp"
// #include "../src/types/IOConfig.hpp"
// #include "../src/types/PrimerConfig.hpp"
// #include "../src/types.hpp"
// #include "../src/utilities.hpp"
//
// // g++ ../PriSeT/apps/barcode_miner_bf.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -lstdc++fs -Wno-unknown-pragmas -lstdc++fs -DNDEBUG -O3 -I~/include -L~/lib -ldivsufsort -o barcode_miner_bf
// // barcode_miner_bf -l <dir_library> -w <dir_work>
//
// using namespace priset;
// using namespace std;
//
// namespace fs = experimental::filesystem;
//
//
// struct Solver
// {
//     vector<size_t> & primerIDs;
//     size_t m;                      // maximal number of primers
//     vector<vector<bool>> & solutions;  // optimal primer combination
//     size_t C_max = 0;              // maximal score
//
//     Solver(vector<size_t> & _primerIDs, int _m, vector<vector<bool>> & _solutions) :
//     primerIDs(_primerIDs), m(_m), solutions(_solutions) {}
//
//     size_t solveBruteForce(vector<bitset<4>> const & C)
//     {
//         solveBruteForceRecursive(C);
//         return C_max;
//     }
//
// private:
//     void solveBruteForceRecursive(vector<bitset<4>> const & C,
//         vector<bool> primerIDs_used = vector<bool>(4, 0), bitset<4> acc = bitset<4>(0), size_t m_tmp = 0)
//     {
//         if (m_tmp < m)
//         {
//             for (size_t primerID : primerIDs)
//             {
//                 if (primerIDs_used[primerID] || ((acc | C.at(primerID)) != acc))
//                     continue;
//                 primerIDs_used[primerID] = 1;
//                 acc |= C.at(primerID);
//                 solveBruteForceRecursive(C, primerIDs_used, acc, ++m_tmp);
//             }
//         }
//         else
//         {
//             if (C_max <= acc.count())
//             {
//                 if (C_max < acc.count()) // delete currently best solutions
//                 {
//                     C_max = acc.count();
//                     solutions.clear();
//                 }
//                 solutions.push_back(primerIDs_used);
//             }
//         }
//     }
// };
//
// int main()
// {
//     int const m = 2; // maximal number of primers
//     vector<size_t> primerIDs = {0, 1, 2, 3};
//     // references
//     unordered_map<size_t, vector<string>> spec2refs;
//     spec2refs[0] = {"AMBERF", "ANBEUF"};
//     spec2refs[1] = {"ANBCODESF", "ANBCODEUF"};
//     spec2refs[2] = {"AMBCPDESF", "ANBEUF"};
//     spec2refs[3] = {"AMBERF", "ANBEUF"};
//
//     // primers
//     vector<pair<string, string>> primers;
//     primers.push_back({"A", "B"});
//     primers.push_back({"C", "D"});
//     primers.push_back({"E", "F"});
//
//     vector<bitset<4>> C;
//
//     for (pair<string, string> primer : primers)
//     {
//         unordered_map<size_t, bitset<4>> amplicon_bitmap;
//         for (auto it = spec2refs.cbegin(); it != spec2refs.cend(); ++it)
//         {
//             auto [sp, refs] = *it;
//             for (string ref : refs)
//             {
//                 auto pos_fwd = ref.find(primer.first);
//                 if (pos_fwd != string::npos)
//                 {
//                     auto pos_rev = ref.find(primer.second);
//                     if (pos_rev != string::npos)
//                     {
//                         string amplicon = ref.substr(pos_fwd+1, pos_rev);
//                         size_t ah = std::hash<string>{}(amplicon);
//                         if (amplicon_bitmap.find(ah) == amplicon_bitmap.end())
//                             amplicon_bitmap[ah] = bitset<4>(0);
//                         amplicon_bitmap[ah].set(sp-1);
//                     }
//                 }
//             }
//         }
//         // all species done, parse amplicon-wise
//         bitset<4> keep(0b1111);
//         bitset<4> c(0);
//         for (auto it = amplicon_bitmap.cbegin(); it != amplicon_bitmap.cend(); ++it)
//         {
//             auto [ah, bv] = *it;
//             size_t set_bits = bv.count();
//             if (set_bits > 1)
//                 keep &= ~bv;
//             else if (set_bits == 1)
//                 c |= bv;
//
//         }
//         // a species can be invalidated afterwards => filter with 'keep'
//         c &= keep;
//         C.push_back(c);
//     }
//
//     // try all combinations
//     vector<vector<bool>> solutions; // primer combination with largest coverage
//     Solver solver = Solver(primerIDs, m, solutions);
//     int score = solver.solveBruteForce(C);
//     cout << "max C: " << score << "\nsolutions: " << endl;
//     for (auto solution : solutions)
//     {
//         cout << "[";
//         for (size_t pID = 1; pID < solution.size(); ++pID)
//         {
//             if (solution[pID])
//                 cout << pID << " ";
//         }
//         cout << "]" << endl;
//     }
//
//
//     return 0;
// }
