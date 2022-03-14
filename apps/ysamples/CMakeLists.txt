add_executable(ysamples  ysamples.cpp)

set_target_properties(ysamples  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(ysamples  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(ysamples PRIVATE yocto)
