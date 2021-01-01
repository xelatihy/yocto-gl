add_executable(yimageproc  yimageproc.cpp)

set_target_properties(yimageproc  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(yimageproc  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yimageproc  yocto)
