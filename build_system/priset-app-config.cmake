# Definitions needed by PriSeT tests and apps independently.

#   PRISET_APP_CXX_FLAGS        -- to be added to CMAKE_CXX_FLAGS

include("${CMAKE_CURRENT_LIST_DIR}/priset-config.cmake")

set (PRISET_APP_CXX_FLAGS "-O3 -DNDEBUG")

# ----------------------------------------------------------------------------
# Macro for printing app messages
# ----------------------------------------------------------------------------

macro(priset_app_message ms)
    message ("${BoldYellow}${msg}${ColourReset}")
endmacro ()
