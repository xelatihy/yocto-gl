add_executable(yconvert  yconvert.cpp)

set_target_properties(yconvert  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(yconvert  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yconvert PRIVATE yocto)
