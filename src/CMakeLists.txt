# g++ priset.cpp -Wno-write-strings -std=c++17 -lstdc++fs -Wall -Wextra -o priset

set(SOURCE_FILES priset.cpp)
set (CMAKE_CXX_FLAGS "-Wno-write-strings -Wall -Wextra")
add_executable(priset ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} stdc++fs)

#message ("${BoldBlue}Linking against genmap ... ${ColourReset}")
#include_directories(../submodules/genmap)

# or:

#ExternalProject_Add (genmap
#    PREFIX "${GENMAP_BUILD_PATH}"
#    DOWNLOAD_COMMAND ""                 # Disable download step, repo was cloned recursively
#    SOURCE_DIR "./submodules/genmap"    # Source dir for configuration
#    INSTALL_DIR "./build/genmap"        # build output directory
#    CMAKE_BUILD_TYPE "Release"
#    GENMAP_NATIVE_BUILD "OFF"
#    CMAKE_CACHE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
#    BUILD_COMMAND ""
#    INSTALL_COMMAND ""
#)

#cmake ./submodules/genmap -DCMAKE_BUILD_TYPE=Release -DGENMAP_NATIVE_BUILD=OFF
#make genmap

#cmake ../submodules/genmap -DGENMAP_NATIVE_BUILD=OFF

#cmake ../submodules/genmap -DCMAKE_BUILD_TYPE=Release -DGENMAP_NATIVE_BUILD=OFF -DEXECUTABLE_OUTPUT_PATH="./bin"
