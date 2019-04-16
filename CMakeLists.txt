# ===========================================================================
#                                 PriSeT
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)

project (priset CXX)

# ----------------------------------------------------------------------------
# Begin of dependency detection
# ----------------------------------------------------------------------------


# Build genmap binary used for FM index building
cmake ../genmap -DCMAKE_BUILD_TYPE=Release -DGENMAP_NATIVE_BUILD=OFF
make genmap
