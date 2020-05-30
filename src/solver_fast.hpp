
#include "solver.hpp"

struct solver_fast : solver
{

    std::vector<std::vector<Result>> results;

    void solve_fast()
    {
        run();

        // set results
        solver::set_results(results);
    }
};
