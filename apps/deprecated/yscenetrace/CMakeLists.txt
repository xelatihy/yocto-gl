add_executable(yscenetrace yscenetrace.cpp)

set_target_properties(yscenetrace PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(yscenetrace PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yscenetrace yocto)
