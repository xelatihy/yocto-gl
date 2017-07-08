cmake_minimum_required (VERSION 3.5)

project (yocto-gl)

option(YOCTO_NO_OPENGL "Skip OpenGL apps" OFF)

add_subdirectory(yocto)
add_subdirectory(apps)
add_subdirectory(tools)
