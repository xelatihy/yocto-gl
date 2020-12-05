add_executable(ymeshtest  ymeshtest.cpp)

set_target_properties(ymeshtest  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(ymeshtest  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(ymeshtest  yocto)

