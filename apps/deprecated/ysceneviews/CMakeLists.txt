add_executable(ysceneviews    ysceneviews.cpp)

set_target_properties(ysceneviews     PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(ysceneviews    PUBLIC ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(ysceneviews    yocto yocto_gui)
