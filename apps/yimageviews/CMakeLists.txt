add_executable(yimageviews   yimageviews.cpp)

set_target_properties(yimageviews    PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(yimageviews   PUBLIC ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yimageviews   yocto yocto_gui)
