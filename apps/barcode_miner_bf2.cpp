// // #include <algorithm>
// // #include <array>
// // #include <chrono>
// // #include <cstdlib>
// // #include <ctime>
// #include <iostream>
// // #include <experimental/filesystem>
// // #include <fstream>
// // #include <numeric>
// // #include <regex>
// // #include <sstream>
// // #include <string>
// // #include <sys/wait.h>
// // #include <unistd.h>
// // #include <vector>
//
// #include "../src/argument_parser.hpp"
// // #include "../src/algorithm.hpp"
// // #include "../src/fm.hpp"
// // #include "../src/IOConfig.hpp"
// // #include "../src/PrimerConfig.hpp"
// #include "../src/solver.hpp"
// // #include "../src/types.hpp"
// // #include "../src/utilities.hpp"
//
// // g++ ../PriSeT/apps/barcode_miner_bf2.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -lstdc++fs -Wno-unknown-pragmas -lstdc++fs -DNDEBUG -O3 -I~/include -L~/lib -ldivsufsort -o barcode_miner_bf2
// // barcode_miner_bf2 -l <dir_library> -w <dir_work> [-m <max_num_primers>] [-s]
//
// using namespace priset;
// using namespace std;
//
// int main(int argc, char ** argv)
// {
//     // set path prefixes for library files
//     IOConfig io_cfg{};
//
//     // Default configuration for primer settings.
//     PrimerConfig primer_cfg{};
//
//     // parse options and init io and primer configurators
//     options opt(argc, argv, primer_cfg, io_cfg);
//
//     // primer_cfg.set_dTm(4);
//     primer_cfg.set_max_primers(2);
//
//
//     // init solver
//     solver_brute_force solver{io_cfg, primer_cfg};
//
//     // solve
//     solver.solve();
//
//     // sort by taxonomic coverage
//     solver.sort_results_by_coverage();
//
//     // output solutions
//     auto [has_next, result] solver.get_next_result():
//     while (has_next)
//     {
//
//         cout << result.to_string() << endl;
//         auto [has_next, result] = solver.get_next_result();
//     }
//     return 0;
// }
