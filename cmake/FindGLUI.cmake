find_path(GLUI_INCLUDE_DIR
    NAMES GL/glui.h
    PATHS /usr/include /usr/local/include
)

find_library(GLUI_LIBRARY
    NAMES glui
    PATHS /usr/lib /usr/local/lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GLUI
    REQUIRED_VARS GLUI_INCLUDE_DIR GLUI_LIBRARY
)

if(GLUI_FOUND)
    set(GLUI_INCLUDE_DIRS ${GLUI_INCLUDE_DIR})
    set(GLUI_LIBRARIES ${GLUI_LIBRARY})
endif()
