# g++ priset.cpp -Wno-write-strings -std=c++17 -lstdc++fs -Wall -Wextra -o priset

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search for zlib as a dependency for SeqAn.
find_package (ZLIB)

# Load the SeqAn module and fail if not found.
find_package (SeqAn REQUIRED)

set (SEQAN_CXX_FLAGS "${SEQAN_CXX_FLAGS} -march=core2")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} -Wall -pedantic -Wno-unknown-pragmas")

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

message ("priset, seqan include_dirs: ${SEQAN_INCLUDE_DIRS}")
message ("priset, SEQAN_LIBRARIES: ${SEQAN_LIBRARIES}")

# PROBLEM: lambda doesn't find calling genmap header

# Add executable and link against SeqAn dependencies.
set(SOURCE_FILES    priset.cpp
                    ../submodules/genmap/src/mappability.hpp
                    ../submodules/genmap/src/common.hpp
                    ../submodules/genmap/src/indexing.hpp
                    ../submodules/genmap/src/mappability.hpp
                    ../submodules/genmap/src/algo.hpp
                    ../submodules/genmap/src/output.hpp
)
# linker error remains with those includes into SOURCE_FILES
#../submodules/genmap/include/lambda/src/mkindex_algo.hpp
#../submodules/genmap/include/lambda/src/mkindex_misc.hpp
#../submodules/genmap/include/lambda/src/mkindex_saca.hpp

add_executable(priset ${SOURCE_FILES})
set_target_properties(priset PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(${PROJECT_NAME} PUBLIC stdc++fs ${SEQAN_LIBRARIES})

#message ("${BoldBlue}Linking against genmap ... ${ColourReset}")
