# ===========================================================================
#                        ▄▄▄·▄▄▄  ▪  .▄▄ · ▄▄▄ .▄▄▄▄▄ 1.0
#                       ▐█ ▄█▀▄ █·██ ▐█ ▀. ▀▄.▀·•██
#                        ██▀·▐▀▀▄ ▐█·▄▀▀▀█▄▐▀▀▪▄ ▐█.▪
#                       ▐█▪·•▐█•█▌▐█▌▐█▄▪▐█▐█▄▄▌ ▐█▌·
#                       .▀   .▀  ▀▀▀▀ ▀▀▀▀  ▀▀▀  ▀▀▀
# ===========================================================================

cmake_minimum_required (VERSION 3.8)

project (priset_app_solver_fast CXX)

set(CMAKE_CXX_STANDARD 17)  # or target_compile_features(priset_lib PUBLIC cxx_std_17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)



# # Add executable and link against SeqAn dependencies.
# set(GLOB_SOURCES ../submodules/genmap/src/*.hpp)
#TODO: define pretty printing macros, for each module its own color
