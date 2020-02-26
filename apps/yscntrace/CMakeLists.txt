add_executable(yscntrace yscntrace.cpp)

set_target_properties(yscntrace PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(yscntrace PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yscntrace yocto)
