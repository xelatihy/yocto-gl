cmake_minimum_required (VERSION 3.12)

project (yocto-gl-ext VERSION 2.0)

add_library(yocto-ext
  stb_image.h stb_image_resize.h stb_image_write.h stb_image.cpp
  cgltf.h cgltf_write.h cgltf.cpp
  filesystem.hpp json.hpp
  tinyexr.h tinyexr.cpp
)

set_target_properties(yocto-ext PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
