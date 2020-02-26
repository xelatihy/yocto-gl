add_executable(yimshproc   yimshproc.cpp)

set_target_properties(yimshproc    PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(yimshproc   PUBLIC ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yimshproc   yocto yocto_gui ${GL_EXTRA_LIBRARIES})
