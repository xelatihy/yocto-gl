cmake_minimum_required (VERSION 3.12)

project (yocto-gl-apps VERSION 2.0)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_executable(yscntrace yscntrace.cpp)
add_executable(yscnproc  yscnproc.cpp)
add_executable(yimgproc  yimgproc.cpp)
add_executable(ymshproc  ymshproc.cpp)

set_target_properties(yscntrace PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
set_target_properties(yscnproc  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
set_target_properties(yimgproc  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
set_target_properties(ymshproc  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

target_link_libraries(yscntrace yocto)
target_link_libraries(yscnproc  yocto)
target_link_libraries(yimgproc  yocto)
target_link_libraries(ymshproc  yocto)

if(YOCTO_EMBREE)
  target_compile_definitions(yscntrace PUBLIC -DYOCTO_EMBREE)
  target_compile_definitions(yscnproc PUBLIC  -DYOCTO_EMBREE)
  target_compile_definitions(yimgproc PUBLIC  -DYOCTO_EMBREE)
  target_compile_definitions(ymshproc PUBLIC  -DYOCTO_EMBREE)
endif(YOCTO_EMBREE)

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

  add_executable(yscnview    yscnview.cpp)
  add_executable(yimshproc   yimshproc.cpp)
  add_executable(yscnitrace  yscnitrace.cpp)
  add_executable(yscnitraces yscnitraces.cpp)
  add_executable(yimgview    yimgview.cpp)
  add_executable(yimgviews   yimgviews.cpp)

  if(YOCTO_EMBREE)
    target_compile_definitions(yscnview PUBLIC -DYOCTO_EMBREE)
    target_compile_definitions(yimshproc PUBLIC  -DYOCTO_EMBREE)
    target_compile_definitions(yscnitrace PUBLIC  -DYOCTO_EMBREE)
    target_compile_definitions(yscnitraces PUBLIC  -DYOCTO_EMBREE)
    target_compile_definitions(yimgview PUBLIC  -DYOCTO_EMBREE)
    target_compile_definitions(yimgviews PUBLIC  -DYOCTO_EMBREE)
  endif(YOCTO_EMBREE)

  set_target_properties(yscnview     PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
  set_target_properties(yimshproc    PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
  set_target_properties(yscnitrace   PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
  set_target_properties(yscnitraces  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
  set_target_properties(yimgview     PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
  set_target_properties(yimgviews    PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

  target_link_libraries(yscnview    yocto yocto-gui ${GL_EXTRA_LIBRARIES})
  target_link_libraries(yimshproc   yocto yocto-gui ${GL_EXTRA_LIBRARIES})
  target_link_libraries(yscnitrace  yocto yocto-gui ${GL_EXTRA_LIBRARIES})
  target_link_libraries(yscnitraces yocto yocto-gui ${GL_EXTRA_LIBRARIES})
  target_link_libraries(yimgview    yocto yocto-gui ${GL_EXTRA_LIBRARIES})
  target_link_libraries(yimgviews   yocto yocto-gui ${GL_EXTRA_LIBRARIES})
endif(YOCTO_OPENGL)
