add_library(yocto
  yocto_math.h yocto_color.h yocto_geometry.h
  yocto_noise.h yocto_sampling.h yocto_shading.h
  yocto_modelio.h yocto_modelio.cpp
  yocto_bvh.h yocto_bvh.cpp
  yocto_shape.h yocto_shape.cpp
  yocto_mesh.h yocto_mesh.cpp 
  yocto_image.h yocto_image.cpp
  yocto_scene.h yocto_scene.cpp
  yocto_trace.h yocto_trace.cpp
  yocto_sceneio.h yocto_sceneio.cpp
  yocto_cli.h
  ext/stb_image.h ext/stb_image_resize.h ext/stb_image_write.h ext/stb_image.cpp
  ext/json.hpp
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
    target_include_directories(yocto PUBLIC "/Program\ Files/Intel/Embree3/include"  "C:/Program\ Files/Intel/Embree3/include")
    target_link_directories(yocto PUBLIC "/Program\ Files/Intel/Embree3/lib" "C:/Program\ Files/Intel/Embree3/lib")
    target_link_libraries(yocto embree3 tbb)
  endif(MSVC)
  if(UNIX AND NOT APPLE)
    target_link_libraries(yocto embree3)
  endif()
endif(YOCTO_EMBREE)

if(YOCTO_DENOISE)
  target_compile_definitions(yocto PUBLIC -DYOCTO_DENOISE)
  if(APPLE)
    target_include_directories(yocto PUBLIC /usr/local/include)
    target_link_libraries(yocto /usr/local/lib/libOpenImageDenoise.dylib)
  endif(APPLE)
  if(MSVC)
    target_include_directories(yocto PUBLIC "/Program\ Files/Intel/Embree3/include"  "C:/Program\ Files/Intel/Embree3/include")
    target_link_directories(yocto PUBLIC "/Program\ Files/Intel/Embree3/lib" "C:/Program\ Files/Intel/Embree3/lib")
    target_link_libraries(yocto embree3 tbb)
  endif(MSVC)
  if(UNIX AND NOT APPLE)
    target_link_libraries(yocto embree3)
  endif()
endif(YOCTO_DENOISE)

# warning flags
if(APPLE)
  target_compile_options(yocto PUBLIC -Wall -Wconversion -Wno-sign-conversion -Wno-implicit-float-conversion)
endif(APPLE)
if(MSVC)
  # target_compile_options(yocto PUBLIC /D_CRT_SECURE_NO_WARNINGS /wd4018 /wd4244 /wd4305 /wd4800 /wd4267)
  target_compile_options(yocto PUBLIC /D_CRT_SECURE_NO_WARNINGS)
  target_compile_options(yocto PUBLIC /EHsc)
endif(MSVC)