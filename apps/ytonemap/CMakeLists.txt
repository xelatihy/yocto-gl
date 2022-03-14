add_executable(ytonemap  ytonemap.cpp)

set_target_properties(ytonemap  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(ytonemap  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(ytonemap PRIVATE yocto)
