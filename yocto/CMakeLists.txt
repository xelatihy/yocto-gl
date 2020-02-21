cmake_minimum_required (VERSION 3.12)

project (yocto-gl-lib VERSION 2.0)

add_library(yocto 
  yocto_math.h yocto_commonio.h 
  yocto_ply.h yocto_obj.h yocto_pbrt.h 
  yocto_shape.h yocto_shape.cpp 
  yocto_image.h yocto_image.cpp
  yocto_trace.h yocto_trace.cpp
  yocto_sceneio.h yocto_sceneio.cpp
  ext/stb_image.h ext/stb_image_resize.h ext/stb_image_write.h ext/stb_image.cpp
  ext/cgltf.h ext/cgltf_write.h ext/cgltf.cpp
  ext/filesystem.hpp ext/json.hpp
  ext/tinyexr.h ext/tinyexr.cpp
)

set_target_properties(yocto PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

if(UNIX AND NOT APPLE)
  find_package(Threads REQUIRED)
  target_link_libraries(yocto Threads::Threads)
endif(UNIX AND NOT APPLE)

if(YOCTO_EMBREE)
  target_compile_definitions(yocto PUBLIC -DYOCTO_EMBREE)
  if(APPLE)
    target_include_directories(yocto PUBLIC /usr/local/include)
    target_link_libraries(yocto /usr/local/lib/libembree3.dylib)
  endif(APPLE)
  if(MSVC)
    target_include_directories(yocto PUBLIC "/Program Files/Intel/Embree3 X64/include" "C:/Program Files/Intel/Embree3 X64/include")
    link_directories("/Program Files/Intel/Embree3 X64/lib" "C:/Program Files/Intel/Embree3 X64/lib")
    target_link_libraries(yocto embree3 tbb tbbmalloc)
  endif(MSVC)
endif(YOCTO_EMBREE)
