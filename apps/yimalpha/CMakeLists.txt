add_executable(yimalpha  yimalpha.cpp)

set_target_properties(yimalpha  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(yimalpha  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yimalpha PRIVATE yocto)
