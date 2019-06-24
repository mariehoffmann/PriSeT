#pragma once

#include <fstream>
#include <regex>
#include <string>
//#include <system>
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

// Compiles Shiny app and starts browser on success.
bool compile_app(priset::io_cfg_type & io_cfg)
{
    char const * cmd = std::string("Rscript " + io_cfg.get_script_file().string() + "\0").c_str();
    //execl("Rscript", &io_cfg.get_script_file().string()[0u], NULL);
    std::string result = exec(cmd);
    std::cout << result << std::endl;
    std::basic_regex const url_rx("Listening on (http\\:\\/\\/[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+\\:[0-9]+)");
    std::smatch url_match;
    if (std::regex_search(result, url_match, url_rx))
    {
        std::cout << "Match found: " << url_match[0].str() << std::endl;
        if (IS_DARWIN)
            execl("open", "-a", "firefox", &url_match[0].str()[0u]); //http://127.0.0.1:3144
        else
            std::cout << "Start app manually in browser: " << url_match[0].str() << std::endl;
    }
    else
    {
        std::cout << "WARNING: Could not extract app URL from " << result << std::endl;
        return false;
    }
    return true;
}

// Generates app script from template.
bool generate_app(priset::io_cfg_type & io_cfg)
{
    // copy first template into newly created `app` folder in working directory
    fs::path const app_template = io_cfg.get_app_template();
    fs::path const script = io_cfg.get_script_file();
    fs::create_directory(script.parent_path());
    std::ifstream src(app_template.string(), std::ios::in);
    std::ofstream dst(script.string(), std::ios::out);

    std::string code((std::istreambuf_iterator<char>(src)), (std::istreambuf_iterator<char>()));
    // replace tags
    std::unordered_map<std::string, std::string> value_map
    {    {"<tax_file>", io_cfg.get_tax_file().string()},
        {"<primer_info>", io_cfg.get_primer_info_file().string()},
        {"<result_file>", io_cfg.get_result_file().string()}
    };
    for (const auto & [tag, value]: value_map)
    {
        code = std::regex_replace(code, std::regex(tag), value);
    }
    dst << code;

    std::cout << "R script copied to: \n" << script << std::endl;
    return true;
}

}  // namespace priset
