add_executable(yimage  yimage.cpp)

set_target_properties(yimage  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(yimage  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yimage  yocto)
