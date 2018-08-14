// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================

/*!\file
 * \brief Core concepts.
 * \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 */

# This file provides functionality common to the different test modules used by
# PriSeT. To build tests, run cmake on one of the sub-folders in this directory
# which contain a CMakeLists.txt.

cmake_minimum_required (VERSION 3.2)

# require PriSeT package
find_package (PriSeT REQUIRED
              HINTS ${CMAKE_CURRENT_LIST_DIR}/../build_system)

include (CheckCXXSourceCompiles)
include (FindPackageHandleStandardArgs)
include (FindPackageMessage)

# ----------------------------------------------------------------------------
# Paths to folders.
# ----------------------------------------------------------------------------

set (SEQAN3_BENCHMARK_CLONE_DIR "${PROJECT_BINARY_DIR}/vendor/benchmark")
set (PRISET_TEST_CLONE_DIR "${PROJECT_BINARY_DIR}/vendor/googletest")

# ----------------------------------------------------------------------------
# Interface targets for the different test modules in priset.
# ----------------------------------------------------------------------------

# priset::test exposes a base set of required flags, includes, definitions and
# libraries which are in common for **all** priset tests
add_library (priset::test INTERFACE IMPORTED)
set_property (TARGET priset::test APPEND PROPERTY INTERFACE_COMPILE_OPTIONS "-pedantic"  "-Wall" "-Wextra" "-Werror")
set_property (TARGET priset::test APPEND PROPERTY INTERFACE_LINK_LIBRARIES "priset::priset" "pthread")
set_property (TARGET priset::test APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${SEQAN3_CLONE_DIR}/test/include/")

# priset::test::performance specifies required flags, includes and libraries
# needed for performance test cases in priset/test/performance
add_library (priset::test::performance INTERFACE IMPORTED)
set_property (TARGET priset::test::performance APPEND PROPERTY INTERFACE_LINK_LIBRARIES "priset::test" "gbenchmark")
set_property (TARGET priset::test::performance APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${SEQAN3_BENCHMARK_CLONE_DIR}/include/")
file(MAKE_DIRECTORY ${SEQAN3_BENCHMARK_CLONE_DIR}/include/) # see cmake bug https://gitlab.kitware.com/cmake/cmake/issues/15052

# priset::test::unit specifies required flags, includes and libraries
# needed for unit test cases in priset/test/unit
add_library (priset::test::unit INTERFACE IMPORTED)
set_property (TARGET priset::test::unit APPEND PROPERTY INTERFACE_LINK_LIBRARIES "priset::test" "gtest_main" "gtest")
set_property (TARGET priset::test::unit APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${PRISET_TEST_CLONE_DIR}/googletest/include/")
file(MAKE_DIRECTORY ${PRISET_TEST_CLONE_DIR}/googletest/include/) # see cmake bug https://gitlab.kitware.com/cmake/cmake/issues/15052

# priset::test::header specifies required flags, includes and libraries
# needed for header test cases in priset/test/header
add_library (priset::test::header INTERFACE IMPORTED)
set_property (TARGET priset::test::header APPEND PROPERTY INTERFACE_LINK_LIBRARIES "priset::test::unit")

# ----------------------------------------------------------------------------
# Commonly shared options for external projects.
# ----------------------------------------------------------------------------

set (SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}")

# ----------------------------------------------------------------------------
# Commonly used macros for the different test modules in priset.
# ----------------------------------------------------------------------------

macro (seqan3_require_ccache)
    find_program (CCACHE_PROGRAM ccache)
    find_package_message (CCACHE_PROGRAM_PRE "Finding program ccache" "[${CCACHE_PROGRAM}]")

    if (NOT CCACHE_PROGRAM)
        find_package_message (CCACHE_PROGRAM "Finding program ccache - Failed" "[${CCACHE_PROGRAM}]")
    else ()
        find_package_message (CCACHE_PROGRAM "Finding program ccache - Success" "[${CCACHE_PROGRAM}]")
        if (CMAKE_VERSION VERSION_LESS 3.4)
            set_property (GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
            set_property (GLOBAL PROPERTY RULE_LAUNCH_LINK "${CCACHE_PROGRAM}")
        else ()
            # New option since cmake >= 3.4:
            # https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_LAUNCHER.html
            set (CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")

            # use ccache in external cmake projects
            list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_COMPILER_LAUNCHER=${CMAKE_CXX_COMPILER_LAUNCHER}")
        endif ()
    endif ()
    unset (CCACHE_PROGRAM)
endmacro ()

macro (seqan3_require_benchmark)
    enable_testing ()

    set (gbenchmark_project_args ${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS})
    list (APPEND gbenchmark_project_args "-DBENCHMARK_ENABLE_TESTING=false")
    # list (APPEND gbenchmark_project_args "-DBENCHMARK_ENABLE_LTO=true")

    include (ExternalProject)
    ExternalProject_Add (
        gbenchmark_project
        PREFIX gbenchmark_project
        GIT_REPOSITORY "https://github.com/google/benchmark.git"
        SOURCE_DIR "${SEQAN3_BENCHMARK_CLONE_DIR}"
        CMAKE_ARGS "${gbenchmark_project_args}"
        UPDATE_DISCONNECTED yes
    )
    unset (gbenchmark_project_args)

    set (gbenchmark_path "${PROJECT_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}benchmark${CMAKE_STATIC_LIBRARY_SUFFIX}")

    add_library (gbenchmark STATIC IMPORTED)
    add_dependencies (gbenchmark gbenchmark_project)
    set_target_properties (gbenchmark PROPERTIES IMPORTED_LOCATION "${gbenchmark_path}")
    set_property (TARGET gbenchmark APPEND PROPERTY INTERFACE_LINK_LIBRARIES "pthread")

    unset (gbenchmark_path)
endmacro ()

macro (seqan3_require_test)
    enable_testing ()

    set (gtest_project_args ${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS})
    list (APPEND gtest_project_args "-DBUILD_GTEST=1")
    list (APPEND gtest_project_args "-DBUILD_GMOCK=0")

    # force that libraries are installed to `lib/`, because GNUInstallDirs might install it into `lib64/`
    list (APPEND gtest_project_args "-DCMAKE_INSTALL_LIBDIR=${PROJECT_BINARY_DIR}/lib/")

    include (ExternalProject)
    ExternalProject_Add (
        gtest_project
        PREFIX gtest_project
        GIT_REPOSITORY "https://github.com/google/googletest.git"
        GIT_TAG "15392f1a38fa0b8c3f13a9732e94b209069efa1c"
        SOURCE_DIR "${PRISET_TEST_CLONE_DIR}"
        CMAKE_ARGS "${gtest_project_args}"
        UPDATE_DISCONNECTED yes
    )
    unset (gtest_project_args)

    # google sets CMAKE_DEBUG_POSTFIX = "d"
    set (gtest_main_path "${PROJECT_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtest_main${CMAKE_STATIC_LIBRARY_SUFFIX}")
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        set (gtest_main_path "${PROJECT_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtest_maind${CMAKE_STATIC_LIBRARY_SUFFIX}")
    endif ()

    add_library (gtest_main STATIC IMPORTED)
    add_dependencies (gtest_main gtest_project)
    set_target_properties (gtest_main PROPERTIES IMPORTED_LOCATION "${gtest_main_path}")

    set (gtest_path "${PROJECT_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtest${CMAKE_STATIC_LIBRARY_SUFFIX}")
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        set (gtest_path "${PROJECT_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtestd${CMAKE_STATIC_LIBRARY_SUFFIX}")
    endif ()

    add_library (gtest STATIC IMPORTED)
    add_dependencies (gtest gtest_main)
    set_target_properties (gtest PROPERTIES IMPORTED_LOCATION "${gtest_path}")
    set_property (TARGET gtest APPEND PROPERTY INTERFACE_LINK_LIBRARIES "pthread")

    unset(gtest_main_path)
    unset(gtest_path)
endmacro ()

macro (add_subdirectories)
    file (GLOB ENTRIES
          RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
          ${CMAKE_CURRENT_SOURCE_DIR}/[!.]*)

    foreach (ENTRY ${ENTRIES})
        if (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${ENTRY})
            if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${ENTRY}/CMakeLists.txt)
                add_subdirectory (${ENTRY})
            endif ()
        endif ()
    endforeach ()
    unset (ENTRIES)
endmacro ()
