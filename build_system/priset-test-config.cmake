# ----------------------------------------------------------------------------
# Macro for printing test messages
# ----------------------------------------------------------------------------
macro(priset_test_message msg)
    message (STATUS "${BoldYellow}${msg}${ColourReset}")
endmacro ()


# ----------------------------------------------------------------------------
# Directories for Google Test Suite and Folder Paths.
# ----------------------------------------------------------------------------
# find_path (PRISET_TEST_WORK_DIR NAMES priset/test/tmp_filename.hpp HINTS "${CMAKE_CURRENT_LIST_DIR}/../test/include/")
# searches for priset-config.cmake
find_package (priset REQUIRED HINTS ${CMAKE_CURRENT_LIST_DIR})

# set and create clone directory for googletest
set (PRISET_TEST_CLONE_DIR "${PROJECT_BINARY_DIR}/vendor/googletest")
# needed for add_library (priset::test::* INTERFACE IMPORTED)
# see cmake bug https://gitlab.kitware.com/cmake/cmake/issues/15052
file(MAKE_DIRECTORY ${PRISET_TEST_CLONE_DIR}/googletest/include/)


include(GoogleTest OPTIONAL)

# Download and build gtest as external project
macro(priset_require_test)
    # populated by cmake call and optional parameter settings
    set (GTEST_PROJECT_CMAKE_ARGS "")
    list (APPEND GTEST_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
    list (APPEND GTEST_PROJECT_CMAKE_ARGS "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
    list (APPEND GTEST_PROJECT_CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}")
    list (APPEND GTEST_PROJECT_CMAKE_ARGS "-DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}")

    # Enables testing for this directory and below.
    enable_testing()
    set (gtest_project_args ${GTEST_PROJECT_CMAKE_ARGS})
    list (APPEND gtest_project_args "-DBUILD_GMOCK=0")

    # force that libraries are installed to `lib/`, because GNUInstallDirs might install it into `lib64/`
    list (APPEND gtest_project_args "-DCMAKE_INSTALL_LIBDIR=${PROJECT_BINARY_DIR}/lib/")

    # google sets CMAKE_DEBUG_POSTFIX = "d"
    set (gtest_main_path "${PROJECT_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtest_main${CMAKE_STATIC_LIBRARY_SUFFIX}")
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        set (gtest_main_path "${PROJECT_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtest_maind${CMAKE_STATIC_LIBRARY_SUFFIX}")
    endif ()

    set (gtest_path "${PROJECT_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtest${CMAKE_STATIC_LIBRARY_SUFFIX}")
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
        set (gtest_path "${PROJECT_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gtestd${CMAKE_STATIC_LIBRARY_SUFFIX}")
    endif ()

    include (ExternalProject)
    ExternalProject_Add (
        gtest_project
        PREFIX gtest_project
        GIT_REPOSITORY "https://github.com/google/googletest.git"
        # we currently have warnings that were introduced in
        # 03867b5389516a0f185af52672cf5472fa0c159c, which are still available
        # in "release-1.8.1", see https://github.com/google/googletest/issues/1419
        GIT_TAG "release-1.10.0"
        SOURCE_DIR "${PRISET_TEST_CLONE_DIR}"
        # INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}"
        CMAKE_ARGS "${gtest_project_args}"
        BUILD_BYPRODUCTS "${gtest_main_path}" "${gtest_path}"
        UPDATE_DISCONNECTED YES
    )

    unset (gtest_project_args)

    add_library (gtest_main STATIC IMPORTED) # IMPORTED: scope of target is this folder and below!
    add_dependencies (gtest_main gtest_project) # ensure gtest_project built before gtest_main
    set_target_properties (gtest_main PROPERTIES IMPORTED_LOCATION "${gtest_main_path}")

    add_library (gtest STATIC IMPORTED)
    add_dependencies (gtest gtest_main)
    set_target_properties (gtest PROPERTIES IMPORTED_LOCATION "${gtest_path}")
    set_property (TARGET gtest APPEND PROPERTY INTERFACE_LINK_LIBRARIES "pthread")

    unset(gtest_main_path)
    unset(gtest_path)
endmacro ()

# Get a specific component of a test file which follows the naming scheme.
# e.g. target_source_file = "types/combine_pattern_test.cpp"
# component:
#  * TARGET_NAME - the target name (e.g. "combine_pattern_test")
#  * TARGET_UNIQUE_NAME - the target name which includes the target_path (e.g. types-combine_pattern_test)
#  * TEST_NAME - the test name which includes the target_path (e.g. "types/combine_pattern")
#  * TARGET_PATH - the path to the target source (e.g. "types")
macro (priset_test_component VAR target_source_file component_name_)
    string (TOUPPER "${component_name_}" component_name)

    get_filename_component (target_relative_path "${target_source_file}" DIRECTORY)
    get_filename_component (target_name "${target_source_file}" NAME_WE)

    if (component_name STREQUAL "TARGET_NAME")
        set (${VAR} "${target_name}")
    elseif (component_name MATCHES "TEST_NAME|TARGET_UNIQUE_NAME")
        if (target_relative_path)
            set (${VAR} "${target_relative_path}/${target_name}")
        else ()
            set (${VAR} "${target_name}")
        endif ()

        if (component_name STREQUAL "TARGET_UNIQUE_NAME")
            string (REPLACE "/" "-" ${VAR} "${${VAR}}")
        endif ()
    elseif (component_name STREQUAL "TARGET_PATH")
        set (${VAR} "${target_relative_path}")
    endif ()

    unset (target_name)
    unset (target_relative_path)
endmacro ()

# loop over all subdirectories to search for tests
macro (add_subdirectories_of directory)
    file (GLOB ENTRIES
          RELATIVE ${directory}
          ${directory}/[!.]*)

    foreach (ENTRY ${ENTRIES})
        if (IS_DIRECTORY ${directory}/${ENTRY})
            if (EXISTS ${directory}/${ENTRY}/CMakeLists.txt)
                add_subdirectory (${directory}/${ENTRY} ${CMAKE_CURRENT_BINARY_DIR}/${ENTRY})
            endif ()
        endif ()
    endforeach ()
    unset (ENTRIES)
endmacro ()


# priset::test exposes a base set of required flags, includes, definitions and
# libraries which are in common to all priset tests
# define interface library for later linking to tests via target_link_libraries()
add_library (priset_test INTERFACE)
target_compile_options (priset_test INTERFACE "-pedantic" "-Wall" "-Wextra" "-Werror" "-Wno-unused-variable")
target_link_libraries (priset_test INTERFACE "priset::priset")
target_include_directories (priset_test INTERFACE "${PRISET_INCLUDE_DIRS}")
add_library (priset::test ALIAS priset_test)

# priset::test::unit specifies required flags, includes and libraries
# needed for unit test cases in priset/test/unit
add_library (priset_test_unit INTERFACE)
# link against Google test suite for unit tests
# set(TEST_UNIT_CXX_FLAG "-lgtest") # obsolete due to next line
target_link_libraries (priset_test_unit INTERFACE "priset::test" "gtest_main" "gtest")
target_include_directories (priset_test_unit INTERFACE "${PRISET_TEST_CLONE_DIR}/googletest/include/")
add_library (priset::test::unit ALIAS priset_test_unit)

macro (add_subdirectories)
    add_subdirectories_of(${CMAKE_CURRENT_SOURCE_DIR})
endmacro ()
