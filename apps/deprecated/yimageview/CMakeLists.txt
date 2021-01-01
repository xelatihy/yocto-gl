add_executable(yimageview    yimageview.cpp)

set_target_properties(yimageview     PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(yimageview    PUBLIC ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yimageview    yocto yocto_gui)
