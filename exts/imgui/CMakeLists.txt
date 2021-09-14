cmake_minimum_required(VERSION 3.12.0)

project(imgui)

find_package(OpenGL REQUIRED)

add_library(imgui STATIC 
  imgui/imgui.cpp
  imgui/imgui_draw.cpp
  imgui/imgui_tables.cpp
  imgui/imgui_widgets.cpp
  imgui/imgui_demo.cpp
  imgui/backends/imgui_impl_glfw.cpp
  imgui/backends/imgui_impl_opengl3.cpp
)

target_include_directories(imgui PUBLIC . imgui imgui/backends)
target_compile_definitions(imgui PRIVATE -DIMGUI_IMPL_OPENGL_LOADER_GLAD)
target_link_libraries(imgui PUBLIC glfw glad ${OPENGL_gl_LIBRARY})

add_library(imgui::imgui ALIAS imgui)

if(UNIX AND NOT APPLE)
target_link_libraries(imgui PUBLIC dl)
endif(UNIX AND NOT APPLE)
