cmake_minimum_required (VERSION 3.5)

project (yocto-gl)

option(YOCTO_NO_OPENGL "Skip OpenGL apps" OFF)

set(CMAKE_CXX_STANDARD 14)
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
if(NOT YOCTO_NO_OPENGL)
    add_definitions(-DYOCTO_OPENGL_APPS)
endif(NOT YOCTO_NO_OPENGL)

add_executable(ysym ysym.cpp yscene.h)
add_executable(ytestgen ytestgen.cpp ext/ArHosekSkyModel.cpp)
add_executable(ytrace ytrace.cpp yscene.h)
add_executable(yobj2gltf yobj2gltf.cpp)
add_executable(ygltfproc ygltfproc.cpp)
add_executable(yimproc yimproc.cpp)

target_link_libraries(ysym yocto)
target_link_libraries(ytestgen yocto)
target_link_libraries(ytrace yocto)
target_link_libraries(yobj2gltf yocto)
target_link_libraries(ygltfproc yocto)
target_link_libraries(yimproc yocto)

if(UNIX AND NOT APPLE)
    set_target_properties(ytrace PROPERTIES COMPILE_FLAGS -pthread LINK_FLAGS -pthread)
endif(UNIX AND NOT APPLE)

if(NOT YOCTO_NO_OPENGL)
    add_executable(yview yview.cpp yscene.h)
    add_executable(yimview yimview.cpp)

    target_link_libraries(yview yocto)
    target_link_libraries(yimview yocto)
endif(NOT YOCTO_NO_OPENGL)
