cmake_minimum_required(VERSION 3.10)

project(Moebius)

set(CMAKE_VERBOSE_MAKEFILE off)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

#set(CMAKE_CXX_FLAGS "-lstdc++")
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
#set(CMAKE_CXX_FLAGS_DEBUG "-g")
#set(CMAKE_CXX_FLAGS_DEBUG "-fsanitize=address")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(CMAKE_CXX_FLAGS "-pthread")
#set(CMAKE_BUILD_TYPE ASAN)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(GLUI REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/inc)
include_directories(${GLUI_INCLUDE_DIRS})

add_library(moebius STATIC
    src/draw.cpp
    src/geometry.cpp
)

target_link_libraries(moebius
    ${GLUI_LIBRARIES}
    OpenGL::GL
    OpenGL::GLU    
    GLUT::GLUT
)

add_subdirectory(examples)
