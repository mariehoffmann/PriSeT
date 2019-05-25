# ===========================================================================
#                                 PriSeT
# ===========================================================================

cmake_minimum_required (VERSION 3.8)
string (ASCII 27 Esc)
set (ColourBold  "${Esc}[1m")
set (ColourReset "${Esc}[m")
set (BoldBlue    "${Esc}[1;34m")

project (priset CXX)

set(CMAKE_CXX_STANDARD 17)  # or target_compile_features(priset_lib PUBLIC cxx_std_17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# This needs to go before "project ()".
#message ("${BoldBlue}Adding submodule genmap ...${ColourReset}")
add_subdirectory(submodules/genmap)
#message ("${BoldBlue}Done.${ColourReset}")

message ("${BoldBlue}Link priset ...${ColourReset}")
add_subdirectory(src)
message ("${BoldBlue}Done.${ColourReset}")
