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

if(YOCTO_OPENGL)
  add_definitions(-DYOCTO_OPENGL)
  find_package(OpenGL REQUIRED)
  set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL " " FORCE)
  set(GLFW_BUILD_TESTS OFF CACHE BOOL " " FORCE)
  set(GLFW_BUILD_DOCS OFF CACHE BOOL " " FORCE)
  set(GLFW_INSTALL OFF CACHE BOOL " " FORCE)
  add_subdirectory(ext/glfw glfw)

  include_directories(.)
  add_definitions(-DIMGUI_IMPL_OPENGL_LOADER_GLAD)

  add_library(yocto_opengl 
    yocto_opengl.h 
    yocto_opengl.cpp 
    ext/imgui/imgui.cpp 
    ext/imgui/imgui_draw.cpp 
    ext/imgui/imgui_widgets.cpp
    ext/imgui/imgui_impl_glfw.cpp 
    ext/imgui/imgui_impl_opengl3.cpp 
    ext/glad/glad.c)
  
  set_target_properties(yocto_opengl PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
  
  target_link_libraries(yocto_opengl glfw ${OPENGL_gl_LIBRARY})

  if(UNIX AND NOT APPLE)
    target_link_libraries(yocto_opengl dl)
  endif(UNIX AND NOT APPLE)

  add_executable(yscnview    yscnview.cpp)
  add_executable(yimshproc   yimshproc.cpp)
  add_executable(yscnitrace  yscnitrace.cpp)
  add_executable(yscnitraces yscnitraces.cpp)
  add_executable(yimgview    yimgview.cpp)
  add_executable(yimgviews   yimgviews.cpp)

  set_target_properties(yscnview     PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
  set_target_properties(yimshproc    PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
  set_target_properties(yscnitrace   PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
  set_target_properties(yscnitraces  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
  set_target_properties(yimgview     PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
  set_target_properties(yimgviews    PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

  target_link_libraries(yscnview    yocto yocto_opengl ${GL_EXTRA_LIBRARIES})
  target_link_libraries(yimshproc   yocto yocto_opengl ${GL_EXTRA_LIBRARIES})
  target_link_libraries(yscnitrace  yocto yocto_opengl ${GL_EXTRA_LIBRARIES})
  target_link_libraries(yscnitraces yocto yocto_opengl ${GL_EXTRA_LIBRARIES})
  target_link_libraries(yimgview    yocto yocto_opengl ${GL_EXTRA_LIBRARIES})
  target_link_libraries(yimgviews   yocto yocto_opengl ${GL_EXTRA_LIBRARIES})
endif(YOCTO_OPENGL)
