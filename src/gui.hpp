#pragma once

#include <fstream>
#include <regex>
#include <string>
#include <system>
#include <unordered_map>

#include "io_cfg_type.hpp"
#include "types.hpp"

#define IS_DARWIN 0

#if defined(__APPLE__) && defined(__MACH__)
    #undef IS_DARWIN
    #define IS_DARWIN 1
#endif

namespace priset::gui
{

// string substitution map to generate valid R-script for frontend
struct TRScriptHelper
{
    using TValueMap = typename std::unordered_map<std::string, std::string>;

    TValueMap value_map{
        {"<tax_file>", io_cfg.get_tax_file()},
        {"<primer_info>", io_cfg.get_primer_info_file()},
        {"<result_file>", io_cfg.get_result_file()}
    };
};

// Compiles Shiny app and starts browser on success.
void compile_app(priset::io_cfg_type & io_cfg)
{
    char const * s = "Rscript " + io_cfg.get_script_file().string() + "\0";
    //execl("Rscript", &io_cfg.get_script_file().string()[0u], NULL);
    std::string result = exec(cmd);
    std::cout << result << std::endl;
    std::basic_regex const url_rx = "Listening on (http\:\/\/\d+\.\d+\.\d+\.\d+:\d+)";
    std::match_results match; // std::smatch matches;
    if (std::regex_search(result.begin(), result.end(), match, url_rx))
    {
        std::cout << "Match found: " << match[0].string() << std::endl;
        if (IS_DARWIN)
            execl("open", "-a", "firefox", &match[0].string()[0u]); //http://127.0.0.1:3144
        else
            std::cout << "Start app manually in browser: " << match[0] << std::endl;
    }
    else
    {
        std::cout << "WARNING: Could not extract app URL from " << result << std::endl;
    }
}

// Generates app script from template.
bool generate_app(priset::io_cfg_type & io_cfg)
{
    TRScriptHelper sh{};

    // copy first template into newly created `app` folder in working directory
    fs::path const template = io_cfg.get_app_template();
    fs::path const script = io_cfg.get_script_file();
    fs::create_directory(fs::path::parent_path(script));
    std::ifstream  src(template.string(), std::ios::in);
    std::ofstream  dst(script.string(), std::ios::out);

    auto buffer = src.rdbuf();
    std::cout << src.rdbuf() << std::endl;
    // replace tags
    for (const auto & [tag, value]: sh.value_map)
    {
        buffer = std::regex_replace(buffer, tag, value);
    }

    dst << buffer;

    std::cout << "R script copied to: \n" << script << std::endl;
    return true;
}

}  // namespace priset
