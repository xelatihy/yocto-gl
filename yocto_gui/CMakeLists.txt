cmake_minimum_required (VERSION 3.12)

project (yocto-gl-gui VERSION 2.0)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(YOCTO_OPENGL)
  add_definitions(-DYOCTO_OPENGL)
  find_package(OpenGL REQUIRED)
  set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL " " FORCE)
  set(GLFW_BUILD_TESTS OFF CACHE BOOL " " FORCE)
  set(GLFW_BUILD_DOCS OFF CACHE BOOL " " FORCE)
  set(GLFW_INSTALL OFF CACHE BOOL " " FORCE)
  add_subdirectory(ext/glfw glfw)

  add_library(yocto-gui 
    yocto_gui.h 
    yocto_gui.cpp 
    ext/imgui/imgui.cpp 
    ext/imgui/imgui_draw.cpp 
    ext/imgui/imgui_widgets.cpp
    ext/imgui/imgui_impl_glfw.cpp 
    ext/imgui/imgui_impl_opengl3.cpp 
    ext/glad/glad.c)
  
  target_include_directories(yocto-gui PRIVATE .)
  target_compile_definitions(yocto-gui PRIVATE -DIMGUI_IMPL_OPENGL_LOADER_GLAD)
  set_target_properties(yocto-gui PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
  
  target_link_libraries(yocto-gui glfw ${OPENGL_gl_LIBRARY})

  if(UNIX AND NOT APPLE)
    target_link_libraries(yocto-gui dl)
  endif(UNIX AND NOT APPLE)
endif(YOCTO_OPENGL)
