add_executable(yscnproc  yscnproc.cpp)

set_target_properties(yscnproc  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(yscnproc  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yscnproc  yocto)
