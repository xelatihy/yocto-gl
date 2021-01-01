add_executable(ysceneview    ysceneview.cpp)

set_target_properties(ysceneview     PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(ysceneview    PUBLIC ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(ysceneview    yocto yocto_gui)
