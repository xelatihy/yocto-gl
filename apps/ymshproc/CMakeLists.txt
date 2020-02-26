add_executable(ymshproc  ymshproc.cpp)

set_target_properties(ymshproc  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(ymshproc  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(ymshproc  yocto)

