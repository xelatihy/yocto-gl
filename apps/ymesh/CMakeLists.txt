add_executable(ymesh  ymesh.cpp)

set_target_properties(ymesh  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(ymesh  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(ymesh  yocto)

if(YOCTO_OPENGL)
target_link_libraries(ymesh  yocto_gui)
endif(YOCTO_OPENGL)
