add_executable(ysculpting    ysculpting.cpp)

set_target_properties(ysculpting     PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(ysculpting    PUBLIC ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(ysculpting    yocto yocto_gui)
