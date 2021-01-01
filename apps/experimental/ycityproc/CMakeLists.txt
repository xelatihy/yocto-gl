add_executable(ycityproc  ycityproc.cpp)

set_target_properties(ycityproc  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(ycityproc  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(ycityproc  yocto)
