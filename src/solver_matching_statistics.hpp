#pragma once

#include "solver.hpp"

struct solver_matchting_statistics : solver
{

    std::vector<std::vector<Result>> results;

    void solve()
    {
        run();

        // set results
        solver::set_results(results);
    }
};
