# searches for priset-config.cmake
find_package (priset REQUIRED HINTS ${CMAKE_CURRENT_LIST_DIR})

# priset::app exposes a base set of required flags, includes, definitions and
# libraries which are in common to all priset apps
# define interface library for later linking to apps via target_link_libraries()
add_library (priset_app INTERFACE)
target_compile_options (priset_app INTERFACE "-pedantic" "-Wall" "-Wextra" "-Werror" "-Wno-unused-variable") # "-ldivsufsort")
target_link_libraries (priset_app INTERFACE "priset::app")
# Link priset_app against more general priset_priset interfaces
target_link_libraries (priset_app INTERFACE priset_priset)

target_include_directories (priset_app INTERFACE "${PRISET_INCLUDE_DIRS}")
add_library (priset::app ALIAS priset_app)

macro(priset_app_macro app_cpp)
    file (RELATIVE_PATH app "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_LIST_DIR}/${app_cpp}")

    get_filename_component(target "${app_cpp}" NAME_WE)
    message("app = ${app}")
    message("target = ${target}")

    add_executable (${target} ${app_cpp})
    target_link_libraries(${target} -lsdsl)
    target_link_libraries(${target} -I$ENV{HOME}/include)
    target_link_libraries(${target} -L$ENV{HOME}/lib)
    
    target_link_libraries (${target} priset::app) 

    unset (app)
    unset (target)
endmacro ()

