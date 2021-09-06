add_executable(yimage-imcli  yimage-imcli.cpp)

set_target_properties(yimage-imcli  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(yimage-imcli  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yimage-imcli  yocto)

if(YOCTO_OPENGL)
target_link_libraries(yimage-imcli  yocto_gui)
endif(YOCTO_OPENGL)
