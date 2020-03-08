add_executable(ysceneitrace  ysceneitrace.cpp)

set_target_properties(ysceneitrace   PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(ysceneitrace  PUBLIC ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(ysceneitrace  yocto yocto_gui)
