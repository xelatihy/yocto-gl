add_library(yocto 
  yocto_math.h yocto_commonio.h 
  yocto_ply.h yocto_obj.h yocto_pbrt.h 
  yocto_bvh.h yocto_bvh.cpp
  yocto_shape.h yocto_shape.cpp 
  yocto_image.h yocto_image.cpp
  yocto_trace.h yocto_trace.cpp
  yocto_sceneio.h yocto_sceneio.cpp
)

set_target_properties(yocto PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

if(UNIX AND NOT APPLE)
  find_package(Threads REQUIRED)
  target_link_libraries(yocto Threads::Threads)
endif(UNIX AND NOT APPLE)

if(YOCTO_EMBREE)
  add_definitions(-DYOCTO_EMBREE)
  if(APPLE)
    include_directories(/usr/local/include)
    link_directories(/usr/local/lib)
    target_link_libraries(yocto libembree3.dylib)
  endif(APPLE)
  if(MSVC)
    include_directories("/Program Files/Intel/Embree3 X64/include" "C:/Program Files/Intel/Embree3 X64/include")
    link_directories("/Program Files/Intel/Embree3 X64/lib" "C:/Program Files/Intel/Embree3 X64/lib")
    target_link_libraries(yocto embree3 tbb tbbmalloc)
  endif(MSVC)
endif(YOCTO_EMBREE)
