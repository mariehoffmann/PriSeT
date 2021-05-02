#
#   PRISET_INCLUDE_DIRS     -- to be passed to include_directories ()
#   PRISET_LIBRARIES        -- to be passed to target_link_libraries ()
#   PRISET_CXX_FLAGS        -- to be added to CMAKE_CXX_FLAGS
#   priset::priset          -- interface target where
#                                  target_link_libraries(target priset::priset)
#                              automatically sets
#                                  target_include_directories(target $SEQAN_INCLUDE_DIRS),
#                                  target_link_libraries(target $SEQAN_LIBRARIES), and
#                                  target_compile_options(target $SEQAN_CXX_FLAGS)
#                              for a target.
#
# ----------------------------------------------------------------------------
# Set paths for PriSeT.
# ----------------------------------------------------------------------------
find_path (PRISET_CLONE_DIR NAMES build_system/priset-test-config.cmake HINTS "${CMAKE_CURRENT_LIST_DIR}/..")
find_path (PRISET_INCLUDE_DIR NAMES version.hpp HINTS "${PRISET_CLONE_DIR}/src")

# ----------------------------------------------------------------------------
# Macros for printing priset and genmap messages
# ----------------------------------------------------------------------------
string (ASCII 27 Esc)
set (ColourBold  "${Esc}[1m")
set (ColourReset "${Esc}[m")
set (BoldGreen   "${Esc}[1;32m")
set (BoldYellow  "${Esc}[1;33m")
set (BoldBlue    "${Esc}[1;34m")
set (BoldMagenta "${Esc}[1;35m")
set (BoldCyan    "${Esc}[1;36m")
set (BoldWhite   "${Esc}[1;37m")

macro(genmap_message msg)
    message ("${BoldBlue}${msg}${ColourReset}")
endmacro ()

macro(priset_status_message msg)
    message (STATUS "${BoldCyan}${msg}${ColourReset}")
endmacro ()

macro(priset_error_message msg)
    message (ERROR "${BoldMagenta}${msg}${ColourReset}")
endmacro ()

macro(priset_fatal_message msg)
    message (FATAL_ERROR "${BoldMagenta}${msg}${ColourReset}")
endmacro ()


# Definitions needed by PriSeT tests and apps independently.
message (STATUS "PRISET_INCLUDE_DIRS = ${PRISET_INCLUDE_DIRS}")

# -Wno-unknown-pragmas to suppress unknown pragma warning in genmap
# -Wno-ignored-qualifiers to suppresses warnings in sdsl,
#   e.g. in int_vector.hpp:1402:8: type qualifiers ignored on function return type
set (PRISET_CXX_FLAGS "-std=c++17 -O3 -lstdc++fs -Wno-unknown-pragmas -Wno-ignored-qualifiers -Wno-deprecated -Wno-error=deprecated-copy -Wno-unused-parameter -lsdsl -ldivsufsort -ldivsufsort64")

# ----------------------------------------------------------------------------
# Find PriSeT include path
# ----------------------------------------------------------------------------

priset_status_message("PRISET_CLONE_DIR = ${PRISET_CLONE_DIR}")
priset_status_message("PRISET_INCLUDE_DIR = ${PRISET_INCLUDE_DIR}")
find_path (PRISET_SUBMODULES_DIR NAMES submodules/genmap HINTS "${PRISET_CLONE_DIR}" "${PRISET_INCLUDE_DIR}/PriSeT")

if (PRISET_INCLUDE_DIR)
    priset_status_message ("PriSeT include dir found:   ${PRISET_INCLUDE_DIR}")
else ()
    priset_fatal_message ("PriSeT include directory could not be found (PRISET_INCLUDE_DIR: '${PRISET_INCLUDE_DIR}')")
endif ()

# ----------------------------------------------------------------------------
# Include check_include_file_cxx macro.
# ----------------------------------------------------------------------------

include (CheckIncludeFileCXX)
include (CheckCXXSourceCompiles)
include (FindPackageHandleStandardArgs)


# # ----------------------------------------------------------------------------
# # Options for CheckCXXSourceCompiles
# # ----------------------------------------------------------------------------
#
# # deactivate messages in check_*
# set (CMAKE_REQUIRED_QUIET       1)
# # use global variables in Check* calls
# set (CMAKE_REQUIRED_INCLUDES    ${CMAKE_INCLUDE_PATH} ${PRISET_INCLUDE_DIR} ${SEQAN3_DEPENDENCY_INCLUDE_DIRS})
# set (CMAKE_REQUIRED_FLAGS       ${CMAKE_CXX_FLAGS})


# ----------------------------------------------------------------------------
# Require C++ Filesystem
# ----------------------------------------------------------------------------

check_include_file_cxx (filesystem PRISET_HAS_FILESYSTEM)
check_include_file_cxx (experimental/filesystem PRISET_HAS_EXP_FILESYSTEM)

if ((NOT PRISET_HAS_FILESYSTEM) AND (NOT PRISET_HAS_EXP_FILESYSTEM))
    priset_fatal_message("PriSeT requires C++17 filesystem support, but the filesystem header was not found.")
endif()


# ----------------------------------------------------------------------------
# Require GenMap and declare interface library.
# ----------------------------------------------------------------------------

add_library (genmap INTERFACE)

target_include_directories(genmap INTERFACE ${PRISET_CLONE_DIR}/submodules/genmap/src)


# ----------------------------------------------------------------------------
# Require SeqAn2 and declare interface library.
# ----------------------------------------------------------------------------

add_library (seqan2 INTERFACE) # declare library interface since genmap is header-only
target_include_directories(seqan2 INTERFACE ${PRISET_CLONE_DIR}/submodules/genmap/include/seqan/include)


# ----------------------------------------------------------------------------
# Require sdsl-lite and declare interface library.
# ----------------------------------------------------------------------------
set(CMAKE_REQUIRED_LIBRARIES "sdsl") 
check_include_file_cxx (sdsl PRISET_HAS_SDSL)

if (PRISET_HAS_SDSL)
    priset_status_message ("Required dependency:        SDSL found.")
else ()
    priset_status_message ("The SDSL library is required, but wasn't found. Get it from https://github.com/xxsds/sdsl-lite")
endif ()


# ----------------------------------------------------------------------------
# Add submodules of dependencies.
# ----------------------------------------------------------------------------

set (PRISET_DEPENDENCY_INCLUDE_DIRS "")

if (PRISET_SUBMODULES_DIR)
    priset_status_message ("PRISET_SUBMODULES_DIR = ${PRISET_SUBMODULES_DIR}")
    file (GLOB submodules ${PRISET_SUBMODULES_DIR}/submodules/*/include)
    foreach (submodule ${submodules})
        if (IS_DIRECTORY ${submodule})
            priset_status_message ("Adding submodule include:  ${submodule}")
            set (PRISET_DEPENDENCY_INCLUDE_DIRS ${PRISET_DEPENDENCY_INCLUDE_DIRS} ${submodule})
        endif ()
    endforeach ()
else ()
    priset_error_message ("Cannot find directory for submodules.")
endif ()


# convert space separated compiler flags into list
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${PRISET_CXX_FLAGS} ${PRISET_LIBRARIES}")
separate_arguments (CMAKE_CXX_FLAGS_LIST UNIX_COMMAND "${CMAKE_CXX_FLAGS}")

add_library (priset_priset INTERFACE)
# target_compile_definitions (priset_priset INTERFACE ${PRISET_DEFINITIONS})
target_compile_options (priset_priset INTERFACE ${CMAKE_CXX_FLAGS_LIST})
target_link_libraries (genmap INTERFACE seqan2)  # necessary?
target_link_libraries (priset_priset INTERFACE genmap)
target_link_libraries (priset_priset INTERFACE seqan2)
# target_link_libraries (priset_priset INTERFACE sdsl)
# target_link_libraries (priset_priset INTERFACE divsufsort)
# target_link_libraries (priset_priset INTERFACE ldivsufsort64)

target_include_directories (priset_priset INTERFACE "${PRISET_INCLUDE_DIR}")
target_include_directories (priset_priset SYSTEM INTERFACE "${PRISET_DEPENDENCY_INCLUDE_DIRS}")
add_library (priset::priset ALIAS priset_priset)

set (PRISET_INCLUDE_DIRS ${PRISET_INCLUDE_DIR} ${PRISET_INCLUDE_DIRS} ${PRISET_DEPENDENCY_INCLUDE_DIRS})
