if(YOCTO_OPENGL)
  add_definitions(-DYOCTO_OPENGL)
  find_package(OpenGL REQUIRED)
  set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL " " FORCE)
  set(GLFW_BUILD_TESTS OFF CACHE BOOL " " FORCE)
  set(GLFW_BUILD_DOCS OFF CACHE BOOL " " FORCE)
  set(GLFW_INSTALL OFF CACHE BOOL " " FORCE)
  add_subdirectory(ext/glfw glfw)

  add_library(ygui 
    yocto_gui.h 
    yocto_gui.cpp 
    ext/imgui/imgui.cpp 
    ext/imgui/imgui_draw.cpp 
    ext/imgui/imgui_widgets.cpp
    ext/imgui/imgui_impl_glfw.cpp 
    ext/imgui/imgui_impl_opengl3.cpp 
    ext/glad/glad.c)
  
  target_include_directories(ygui PRIVATE .)
  target_include_directories(ygui PUBLIC ${CMAKE_SOURCE_DIR}/lib)
  target_compile_definitions(ygui PRIVATE -DIMGUI_IMPL_OPENGL_LOADER_GLAD)
  set_target_properties(ygui PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
  
  target_link_libraries(ygui glfw ${OPENGL_gl_LIBRARY})

  if(UNIX AND NOT APPLE)
    target_link_libraries(ygui dl)
  endif(UNIX AND NOT APPLE)
endif(YOCTO_OPENGL)
