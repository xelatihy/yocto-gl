add_executable(yshade  yshade.cpp)

set_target_properties(yshade  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(yshade  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yshade  yocto)

if(YOCTO_OPENGL)
target_link_libraries(yshade  yocto_gui)
endif(YOCTO_OPENGL)
