add_executable(yimshproc   yimshproc.cpp)

set_target_properties(yimshproc    PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(yimshproc   PUBLIC ${CMAKE_SOURCE_DIR}/lib)
target_link_libraries(yimshproc   ygl ygui ${GL_EXTRA_LIBRARIES})
