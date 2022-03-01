add_executable(ycutrace  ycutrace.cpp)

set_target_properties(ycutrace  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(ycutrace  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(ycutrace PRIVATE yocto)
