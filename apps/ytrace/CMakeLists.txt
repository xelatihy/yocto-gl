add_executable(ytrace  ytrace.cpp)

set_target_properties(ytrace  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(ytrace  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(ytrace  yocto)

if(YOCTO_OPENGL)
target_link_libraries(ytrace  yocto_gui)
endif(YOCTO_OPENGL)
