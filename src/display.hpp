#pragma once

#include <fstream>
#include <string>

#include "io_config.hpp"
#include "types.hpp"

namespace priset
{

// call Shiny to display graphically primer candidates

void display(priset::io_config & io_cfg)
{
    // write R script based on template
    fs::path script_path = io_cfg.get_work_dir() / "display";
    // create directory
    //mkdir
    std::ifstream file(script_path);
    std::string str;
    while (std::getline(file, str))
    {
        // Process str
        // TODO: continue here, write back modified lines
    }
    // Run R script
    // Rscript ../PriSeT/gui/app.R
}

}  // namespace priset
