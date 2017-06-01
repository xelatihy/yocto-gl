cmake_minimum_required (VERSION 3.5)

project (yocto-gl)

option(INCLUDE_OPENGL_DEPS "Include OpenGL dependencies and apps" OFF)

add_subdirectory(yocto)
add_subdirectory(apps)
