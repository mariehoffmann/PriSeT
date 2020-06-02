#include <iostream>

#include "../src/argument_parser.hpp"
#include "../src/solver.hpp"
#include "../src/types/all.hpp"

// g++ ../PriSeT/apps/plankton_denovo2.cpp -Wno-write-strings -std=c++17 -Wall -Wextra -lstdc++fs -Wno-unknown-pragmas -lstdc++fs -DNDEBUG -O3 -I~/include -L~/lib -ldivsufsort -o plankton_denovo2
// plankton_denovo2 -l <dir_library> -w <dir_work> [-m <max_num_primers>] [-s]

using namespace priset;
using namespace std;

int main(int argc, char ** argv)
{
    // set path prefixes for library files
    IOConfig io_cfg{};

    // Default configuration for primer settings.
    PrimerConfig primer_cfg{};

    // parse options and init io and primer configurators
    options opt(argc, argv, primer_cfg, io_cfg);

    // primer_cfg.set_dTm(4);
    primer_cfg.set_dTm(2);

    // init solver
    solver_fast solver{io_cfg, primer_cfg};

    // solve
    solver.solve();

    // sort results by frequency.
    solver.sort_results_by_frequency();

    // output solutions
    cout << solver.get_header() << endl;
    std::vector<Result> results;
    auto [has_next, results] solver.get_next_result():
    size_t ctr{0};
    while (has_next)
    {
        cout << "Result #" << ++ctr << ":\n";
        for (Result result : results)
            cout << "\t" << result.to_string() << endl;
        auto [has_next, results] = solver.get_next_result();
    }
    return 0;
}
