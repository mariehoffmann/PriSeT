#pragma once

#include <vector>

#include "Solver.hpp"
#include "types/Result.hpp"

namespace priset
{

// struct solver;

struct SolverFast : Solver
{
    // solver_fast() = default;

    void solve()
    {
        Solver::run();
        // wrap single pairs as groups
        if (Solver::primer_cfg.get_primer_set_size() == 1)
            Solver::as_groups();
        // greedy: group primer pairs with largest coverage
        else
            Solver::group_by_max_coverage();

    }
};

}  // namespace priset
