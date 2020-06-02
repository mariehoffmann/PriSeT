
#include "solver.hpp"

struct solver_fast : solver
{

    std::vector<std::vector<Result>> results;
    void solve_fast()
    {
        run();
        // wrap single pairs as groups
        if (solver::primer_cfg.get_primer_set_size() == 1)
            solver::as_groups();
    }
};
