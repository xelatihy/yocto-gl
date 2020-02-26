add_executable(yimgview    yimgview.cpp)

set_target_properties(yimgview     PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(yimgview    PUBLIC ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yimgview    yocto yocto_gui)
