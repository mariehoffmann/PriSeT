
#include "solver.hpp"

struct solver_fast : solver
{

    std::vector<std::vector<TResult>> results;

    void solve_fast()
    {
        run();

        // set results
        solver::set_results(results);
    }
};
