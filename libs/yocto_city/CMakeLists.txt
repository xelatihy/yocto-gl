add_executable(yocto_city  yocto_city.cpp)

set_target_properties(yocto_city  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(yocto_city  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yocto_city  yocto yocto_gui)
