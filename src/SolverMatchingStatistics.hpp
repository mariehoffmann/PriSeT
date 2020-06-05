#pragma once

#include "Solver.hpp"

struct SolverMatchtingStatistics : Solver
{

    std::vector<std::vector<Result>> results;

    void solve()
    {
        run();

        // set results
        Solver::set_results(results);
    }
};
