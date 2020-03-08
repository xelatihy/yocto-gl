add_executable(yshapeproc  yshapeproc.cpp)

set_target_properties(yshapeproc  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(yshapeproc  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yshapeproc  yocto)

