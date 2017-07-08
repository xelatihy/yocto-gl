cmake_minimum_required (VERSION 3.5)

project (yocto-gl)

option(YOCTO_NO_OPENGL "Skip OpenGL apps" OFF)
option(YOCTO_CPP11 "Compile as C++11" OFF)

if(NOT YOCTO_CPP11)
set(CMAKE_CXX_STANDARD 14)
else(NOT YOCTO_CPP11)
set(CMAKE_CXX_STANDARD 11)
endif(NOT YOCTO_CPP11)

set(CMAKE_CXX_STANDARD_REQUIRED on)
# set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# compile flags
if(APPLE)
    add_definitions(-Wno-missing-braces)
endif(APPLE)
if(WIN32)
    add_definitions(/D_CRT_SECURE_NO_WARNINGS /wd4018 /wd4244 /wd4305 /wd4800 /wd4267)
endif(WIN32)

add_executable(ygltfgen ygltfgen.cpp ../yocto/yocto_utils.h ../yocto/yocto_utils.cpp)
add_executable(ydocgen ydocgen.cpp ../yocto/yocto_utils.h ../yocto/yocto_utils.cpp)
