#pragma once

#include <fstream>
#include <regex>
#include <string>
#include <system>
#include <unordered_map>

#include "io_cfg_type.hpp"
#include "types.hpp"

namespace priset
{

template<typename TAppValueMap>
void fill_value_map(TAppValueMap & value_map)
{
    // TODO: more tags here
    value_map["<tax_file>"] = io_cfg.get_tax_file();

}

// call Shiny to display graphically primer candidates

void generate_app(priset::io_cfg_type & io_cfg)
{
    // fill value map
    using TAppValueMap = typename std::unordered_map<std::string, std::string>;
    TAppValueMap app_values;
    fill_value_map<TAppValueMap>(value_map);

    // path to generated template
    fs::path script_file = io_cfg.get_work_dir() / "app" / "app.R";

    ostream outFile(script_file);
    istream readFile(io_cfg.get_app_template());

    // create directory
    //mkdir
    std::ifstream file(script_path);
    std::string line;
    std::string buffer = "";
    while (std::getline(file, line))
    {
        buffer += line;
    }
    // replace tags
    for (const auto & [tag, value]: value_map)
    {

        buffer = std::regex_replace(buffer, tag, value);

    }
    std::cout << "new file: \n" << buffer << std::endl;
//    outFile << buffer;
    // Run R script
    // Rscript ../PriSeT/gui/app.R
    //system("Rscript {}")
}

}  // namespace priset
