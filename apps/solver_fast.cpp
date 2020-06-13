#include <iostream>

#include "../src/argument_parser.hpp"
#include "../src/SolverFast.hpp"
#include "../src/types/all.hpp"

/*
 * Compile this app with
 *  g++ ../PriSeT/apps/solver_fast.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -lstdc++fs -Wno-unknown-pragmas -DNDEBUG -O3 -I ~/include -L ~/lib -lsdsl -ldivsufsort -ldivsufsort64 -g -o solver_fast
* and execute with
*   solver_fast -l <dir_library> -w <dir_work> [-m <max_num_primers>] [-s]
*/

int main(int argc, char ** argv)
{
    // set path prefixes for library files
    priset::IOConfig io_cfg{};

    // Default configuration for primer settings.
    priset::PrimerConfig primer_cfg{};

    // parse options and init io and primer configurators
    priset::options opt(argc, argv, primer_cfg, io_cfg);

    // Set melting temperature difference constraint other than default (5 Kelvin).
    primer_cfg.set_primer_dTm(4);

    // init solver
    priset::SolverFast solver{io_cfg, primer_cfg};

    // // solve
    // solver.solve();
    //
    // // sort results by frequency.
    // solver.sort_results_by_frequency();
    //
    // // output solutions
    // cout << solver.get_header() << endl;
    // std::vector<Result> results;
    // auto [has_next, results] solver.get_next_result():
    // size_t ctr{0};
    // while (has_next)
    // {
    //     cout << "Result #" << ++ctr << ":\n";
    //     for (Result result : results)
    //         cout << "\t" << result.to_string() << endl;
    //     auto [has_next, results] = solver.get_next_result();
    // }
    return 0;
}
