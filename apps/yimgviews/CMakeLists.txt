add_executable(yimgviews   yimgviews.cpp)

set_target_properties(yimgviews    PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(yimgviews   PUBLIC ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yimgviews   yocto yocto_gui)
