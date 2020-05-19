# run with cmake . &&
cmake_minimum_required(VERSION 2.8.11)

# https://foonathan.net/2016/07/cmake-dependency-handling/

set(CMAKE_CXX_COMPILER "/opt/local/bin/g++")

project(UnitTests)

if(POLICY CMP0022)
  cmake_policy(SET CMP0022 NEW)
endif()

if(POLICY CMP0054)
  cmake_policy(SET CMP0054 NEW)
endif()


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

################################## GTest ######################################
# https://github.com/google/googletest/blob/master/googletest/README.md
# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
  include_directories("${gtest_SOURCE_DIR}/include")
endif()


################################## SeqAn ######################################
# Find SeqAn with hint for search.
set (SEQAN_FIND_DEPENDENCIES ZLIB)
find_package (SeqAn REQUIRED)

# Add include directories.
include_directories(${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

##################################### SDSL #####################################
#include_directories(/Users/troja/include)
message(WARNING ${PROJECT_SOURCE_DIR})
# include_directories(${PROJECT_SOURCE_DIR}/submodules/sdsl-lite/)
include_directories(${PROJECT_SOURCE_DIR}/../submodules/sdsl-lite/include/sdsl)

# add_subdirectory(${PROJECT_SOURCE_DIR}/../submodules/sdsl-lite/include)
# -isystem /Users/troja/git/PriSeT/submodules/sdsl-lite/include
add_compile_options("-Icore/include -Iextras/include -I/opt/local/include -pedantic -Wno-write-strings -W -Wall -Wextra -DNDEBUG -O3 -L/Users/troja/lib -ldivsufsort -ldivsufsort64 -lstdc++fs -lrt -lstdc++fs -Wno-dev -Wno-long-long -Wno-variadic-macros -DNDEBUG -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0")


################################# Add Tests ####################################
# include(GoogleTest)
# Now simply link against gtest or gtest_main as needed. Eg
add_executable(dna_test dna_test.cpp)
target_link_libraries(dna_test gtest_main ${SEQAN_LIBRARIES} ${GTEST_LIBRARY})
add_test(NAME dna_test COMMAND dna_test)

#set_tests_properties(${types_tests}   PROPERTIES TIMEOUT 10)
#
# target_link_libraries(dna_test
#     PUBLIC ${SEQAN_LIBRARIES}
#     ${GTEST_LIBRARY}
# )
