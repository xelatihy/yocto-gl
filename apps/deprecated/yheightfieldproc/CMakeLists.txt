add_executable(yheightfieldproc  yheightfieldproc.cpp)

set_target_properties(yheightfieldproc  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(yheightfieldproc  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yheightfieldproc  yocto)

