


struct SolverBruteForce : Solver
{
private:
    vector<bitset<4>> C,
    vector<bool> primerIDs_used = vector<bool>(4, 0);
    bitset<4> acc = bitset<4>(0);
    size_t m_tmp = 0;





    void solve_brute_force_recursive()
    {

        if (m_tmp < m)
        {
            for (size_t primerID : primerIDs)
            {
                if (primerIDs_used[primerID] || ((acc | C.at(primerID)) != acc))
                    continue;
                primerIDs_used[primerID] = 1;
                acc |= C.at(primerID);
                solveBruteForceRecursive(C, primerIDs_used, acc, ++m_tmp);
            }
        }
        else
        {
            if (C_max <= acc.count())
            {
                if (C_max < acc.count()) // delete currently best solutions
                {
                    C_max = acc.count();
                    solutions.clear();
                }
                solutions.push_back(primerIDs_used);
            }
        }
    }

public:
    void solve_brute_force()
    {
        run();
        compute_C():
        solve_brute_force_recursive();
    }

};
