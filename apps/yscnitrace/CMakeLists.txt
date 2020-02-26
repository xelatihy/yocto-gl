add_executable(yscnitrace  yscnitrace.cpp)

set_target_properties(yscnitrace   PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(yscnitrace  PUBLIC ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yscnitrace  yocto yocto_gui)
