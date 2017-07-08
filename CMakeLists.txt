cmake_minimum_required (VERSION 3.5)

project (yocto-gl)

option(YOCTO_NO_OPENGL "Skip OpenGL apps" OFF)
option(YOCTO_CPP11 "Compile as C++11" OFF)

add_subdirectory(yocto)
add_subdirectory(apps)
add_subdirectory(tools)
