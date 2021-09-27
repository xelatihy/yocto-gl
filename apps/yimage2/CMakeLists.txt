add_executable(yimage2  yimage2.cpp)

set_target_properties(yimage2 PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(yimage2  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yimage2 PRIVATE yocto)
