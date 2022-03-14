add_executable(yconverts  yconverts.cpp)

set_target_properties(yconverts  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(yconverts  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yconverts PRIVATE yocto)
