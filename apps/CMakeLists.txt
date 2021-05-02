cmake_minimum_required(VERSION 3.7)
project(priset_app CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# searches for priset-app-config.cmake, which will also include priset-config.cmake
find_package (priset-app REQUIRED HINTS ${CMAKE_CURRENT_LIST_DIR}/../build_system)

include (CheckCXXSourceCompiles)
include (FindPackageHandleStandardArgs)
include (FindPackageMessage)

priset_app_macro(solver_fast.cpp)
