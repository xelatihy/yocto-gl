add_executable(yscnview    yscnview.cpp)

set_target_properties(yscnview     PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(yscnview    PUBLIC ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yscnview    yocto yocto_gui)
