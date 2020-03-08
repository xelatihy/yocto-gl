add_executable(ysceneproc  ysceneproc.cpp)

set_target_properties(ysceneproc  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(ysceneproc  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(ysceneproc  yocto)
