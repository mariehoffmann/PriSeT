// ============================================================================
//                    PriSeT - The Primer Search Tool
// ============================================================================

cmake_minimum_required (VERSION 3.2)
project (priset_test_unit CXX)

include (../priset-test.cmake)

macro (priset_test target_cpp)
    # $target_cpp = "pod_tuple_test.cpp"
    #   * will register the global TARGET name "pod_tuple_test" and
    #   * will register the test case name "core_pod_tuple" if
    #     pod_tuple_test.cpp is in test/core/
    #
    # NOTE(marehr): ".+/test/" REGEX is greedy, that means
    # /test/test/test/hello_test.cpp will result in an empty `test_path`
    string (REGEX REPLACE "_test.cpp$" "" target_name ${target_cpp})
    string (REGEX REPLACE ".+/test/" "" test_path ${CMAKE_CURRENT_LIST_DIR})
    set (target "${target_name}_test")

    add_executable (${target} ${target_cpp})
    target_link_libraries (${target} priset::test::unit)
    add_test (NAME "${test_path}/${target}" COMMAND ${target})

    unset (target_name)
    unset (test_path)
    unset (target)
endmacro ()

priset_require_ccache ()
priset_require_test ()

add_subdirectories ()
