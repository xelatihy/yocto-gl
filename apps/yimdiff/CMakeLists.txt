add_executable(yimdiff  yimdiff.cpp)

set_target_properties(yimdiff  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(yimdiff  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yimdiff PRIVATE yocto)
