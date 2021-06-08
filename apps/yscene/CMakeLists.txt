add_executable(yscene  yscene.cpp)

set_target_properties(yscene  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(yscene  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yscene  yocto)

if(YOCTO_OPENGL)
target_link_libraries(yscene  yocto_gui)
endif(YOCTO_OPENGL)
