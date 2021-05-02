#pragma once

#include "Solver.hpp"
#include "types/IOConfig.hpp"
#include "types/PrimerConfig.hpp"
#include "types/Result.hpp"

namespace priset
{

struct SolverFast : Solver
{
    SolverFast(IOConfig & _io_cfg, PrimerConfig & _primer_cfg) :
        Solver(_io_cfg, _primer_cfg) {}

    void run()
    {
        Solver::run();
        // wrap single pairs as groups
        if (Solver::primer_cfg.get_primer_set_size() == 1)
        {
            Solver::as_groups();
            // Solver::sort_groups_by_coverage();
            Solver::sort_results_by_frequency();
        }
        // greedy: group primer pairs with largest coverage
        else
            Solver::group_by_max_coverage_exact();

    }
};

}  // namespace priset
