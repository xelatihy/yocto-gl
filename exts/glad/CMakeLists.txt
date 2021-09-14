cmake_minimum_required(VERSION 3.12.0)

project(glad)

find_package(OpenGL REQUIRED)

add_library(glad STATIC 
  glad/glad.c
)

target_include_directories(glad PUBLIC . glad)
target_link_libraries(glad PUBLIC ${OPENGL_gl_LIBRARY})

add_library(glad::glad ALIAS glad)

if(UNIX AND NOT APPLE)
target_link_libraries(glad PUBLIC dl)
endif(UNIX AND NOT APPLE)
