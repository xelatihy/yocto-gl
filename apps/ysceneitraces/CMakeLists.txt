add_executable(ysceneitraces ysceneitraces.cpp)

set_target_properties(ysceneitraces  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(ysceneitraces PUBLIC ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(ysceneitraces yocto yocto_gui)
