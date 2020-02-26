add_executable(yimgproc  yimgproc.cpp)

set_target_properties(yimgproc  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(yimgproc  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yimgproc  yocto)
