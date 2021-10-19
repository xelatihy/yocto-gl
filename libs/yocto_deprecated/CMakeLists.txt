add_library(yocto_deprecated STATIC)
target_sources(yocto_deprecated PRIVATE 
  yocto_json.h yocto_json.cpp 
  yocto_cli.h yocto_cli.cpp 
  ext/json.hpp
)

set_target_properties(yocto_deprecated PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

target_include_directories(yocto_deprecated PRIVATE ext/)
target_include_directories(yocto_deprecated PUBLIC ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yocto_deprecated PUBLIC yocto)

if(YOCTO_OPENGL)
  find_package(OpenGL REQUIRED)

  target_sources(yocto_deprecated PRIVATE 
    yocto_opengl.h 
    yocto_opengl.cpp 
    yocto_shade.h 
    yocto_shade.cpp 
    yocto_imgui.h 
    yocto_imgui.cpp
    yocto_glview.h 
    yocto_glview.cpp)
  
  target_include_directories(yocto_deprecated PRIVATE .)
  target_include_directories(yocto_deprecated PRIVATE ./ext)
  target_include_directories(yocto_deprecated PUBLIC ${CMAKE_SOURCE_DIR}/libs)
  target_compile_definitions(yocto_deprecated PRIVATE -DIMGUI_IMPL_OPENGL_LOADER_GLAD)
  target_compile_definitions(yocto_deprecated PUBLIC -DYOCTO_OPENGL)
  set_target_properties(yocto_deprecated PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

  target_link_libraries(yocto_deprecated PUBLIC yocto)

  target_include_directories(yocto_deprecated PRIVATE ${CMAKE_SOURCE_DIR}/exts/glfw/include)
  target_link_libraries(yocto_deprecated PUBLIC glad imgui glfw ${OPENGL_gl_LIBRARY})

  if(UNIX AND NOT APPLE)
    target_link_libraries(yocto_deprecated PUBLIC dl)
  endif(UNIX AND NOT APPLE)
endif(YOCTO_OPENGL)

# warning flags
if(APPLE)
  target_compile_options(yocto_deprecated PUBLIC -Wall -Wconversion -Wno-sign-conversion -Wno-implicit-float-conversion)
endif(APPLE)
if(MSVC)
  target_compile_options(yocto_deprecated PUBLIC /D_CRT_SECURE_NO_WARNINGS /wd4018 /wd4244 /wd4305 /wd4800 /wd4267)
  target_compile_options(yocto_deprecated PUBLIC /EHsc)
endif(MSVC)
