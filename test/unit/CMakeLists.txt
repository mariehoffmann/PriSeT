
cmake_minimum_required (VERSION 3.7)
project(priset_test_unit CXX)

# ----------------------------------------------------------------------------
# Required packages
# ----------------------------------------------------------------------------
# searches for priset-config.cmake
find_package (priset-test REQUIRED HINTS ${CMAKE_CURRENT_LIST_DIR}/../../build_system)
# include (${CMAKE_CURRENT_LIST_DIR}/../../build_system/priset-test-config.cmake)

include (CheckCXXSourceCompiles)
include (FindPackageHandleStandardArgs)
include (FindPackageMessage)

# priset_require_ccache ()
priset_require_test()


# macro for adding and linking unit tests
macro(priset_unit_test unit_test_cpp)
    file (RELATIVE_PATH unit_test "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_LIST_DIR}/${unit_test_cpp}")

    priset_test_component (target "${unit_test}" TARGET_NAME)
    priset_test_component (test_name "${unit_test}" TEST_NAME)

    message("unit_test = ${unit_test}")
    message("target = ${target}")
    message("test_name = ${test_name}")

    add_executable (${target} ${unit_test_cpp})
    target_link_libraries (${target} priset::test::unit)  # TODO: check which headers are in this namespace
    if (NOT CMAKE_VERSION VERSION_LESS 3.10) # cmake >= 3.10
        gtest_discover_tests(${target} TEST_PREFIX "${test_name}::" PROPERTIES TIMEOUT "30")
    else ()
        add_test (NAME "${test_name}" COMMAND ${target})
    endif ()

    unset (unit_test)
    unset (target)
    unset (test_name)
endmacro ()

add_subdirectory(filter)
add_subdirectory(types)
