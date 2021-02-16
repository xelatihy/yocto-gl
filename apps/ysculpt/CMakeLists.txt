add_executable(ysculpt  ysculpt.cpp)

set_target_properties(ysculpt  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(ysculpt  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(ysculpt  yocto)

if(YOCTO_OPENGL)
target_link_libraries(ysculpt  yocto_gui)
endif(YOCTO_OPENGL)
