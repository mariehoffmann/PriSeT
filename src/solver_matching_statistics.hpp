#pragma once

#include "solver.hpp"

struct solver_matchting_statistics : solver
{

    std::vector<std::vector<TResult>> results;

    void solve()
    {
        run();

        // set results
        solver::set_results(results);
    }
};
