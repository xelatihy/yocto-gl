add_executable(yscnitraces yscnitraces.cpp)

set_target_properties(yscnitraces  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(yscnitraces PUBLIC ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yscnitraces yocto yocto_gui)
