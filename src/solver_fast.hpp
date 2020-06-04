#pragma once

#include <vector>

#include "solver.hpp"
#include "types/Result.hpp"

namespace priset
{

// struct solver;

struct solver_fast : solver
{
    solver_fast() = default;

    void solve()
    {
        solver::run();
        // wrap single pairs as groups
        if (solver::primer_cfg.get_primer_set_size() == 1)
            solver::as_groups();
        // greedy: group primer pairs with largest coverage
        else
            solver::group_by_max_coverage();
        
    }
};

}  // namespace priset
