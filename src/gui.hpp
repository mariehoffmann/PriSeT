#pragma once

#include <fstream>
#include <regex>
#include <string>
#include  <sys/types.h>
//#include <system>
#include <unistd.h>
#include <unordered_map>


#include "io_cfg_type.hpp"
#include "types.hpp"
#include "utilities.hpp"

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
    std::cout << "Enter compile app\n";
    std::string cmd_str = "Rscript " + io_cfg.get_script_file().string() + "\0";
    std::cout << "cmd_str = " << cmd_str << std::endl;
    char const * cmd = cmd_str.c_str();
    std::cout << "cmd = " << std::string(cmd) << std::endl;
    //exit(0);
    pid_t pid;

    switch(pid = fork()) {
        case -1:   /* error for fork() */
            break;
        case 0:   /* child process */
            {
                std::cout << "Compile script ...\n";
                exec(cmd);
                break;
            }
        default:   /* parent process */
            {
                sleep(2);
                //std::cout << "Launch script ...\n";
                //pid_t pid2 = fork();
                //if (!pid2) exec(cmd2);
                break;
            }
    }
    return true;
}

// Generates app script from template.
bool generate_app(priset::io_cfg_type & io_cfg)
{
    // copy first template into newly created `app` folder in working directory
    fs::path const app_template = io_cfg.get_app_template();
    std::cout << "app_template: " << app_template << std::endl;
    assert(fs::exists(app_template) == true);
    fs::path const script = io_cfg.get_script_file();
    std::cout << "script path = " << script << std::endl;
    fs::create_directory(script.parent_path());
    std::ifstream src(app_template.string(), std::ios::in);
    std::ofstream dst(script.string(), std::ios::out);

    std::string code;

    src.seekg(0, std::ios::end);
    code.reserve(src.tellg());
    src.seekg(0, std::ios::beg);

    code.assign((std::istreambuf_iterator<char>(src)), std::istreambuf_iterator<char>());

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

    //std::cout << "code = " << code << std::endl;
    std::cout << "R script copied to: \n" << script << std::endl;
    return true;
}

}  // namespace priset
