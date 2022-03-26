add_executable(yshape  yshape.cpp)

set_target_properties(yshape  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(yshape  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yshape PRIVATE yocto)
