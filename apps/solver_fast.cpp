#include <iostream>

#include "../src/argument_parser.hpp"
#include "../src/SolverFast.hpp"
#include "../src/types/all.hpp"

/*
*   Example application that creates an index on a given FASTA file by running first 
*   the compiled binary:
*       ./solver_fast -i -l <path_to_fasta> -w <work_dir|idx_dir>
*   and then computes primer pairs:
*       ./solver_fast -s -e exp.cfg -w <idx_dir>
*/
int main(int argc, char ** argv)
{
    // set path prefixes for library files
    priset::IOConfig io_cfg{};

    // Default configuration for primer settings.
    priset::PrimerConfig primer_cfg{};

    // parse options and init io and primer configurators
    priset::options opt(argc, argv, primer_cfg, io_cfg);

    // Set melting temperature difference constraint other than default (5 Kelvin).
    // primer_cfg.set_primer_dTm(4);
    // primer_cfg.set_primer_min_len(16);
    // primer_cfg.set_primer_max_len(25);

    // init solver
    priset::SolverFast solver{io_cfg, primer_cfg};

    // run solver
    solver.run();
    
    // return if only task was to compute index
    if (!io_cfg.get_skip_idx())
        return 0;

    // sort results by frequency.
    solver.sort_results_by_frequency();
    
    // output solutions
    std::cout << solver.get_header() << std::endl;
    std::string s = solver.generate_table();
    std::cout << s << std::endl;
    
    return 0;
}
