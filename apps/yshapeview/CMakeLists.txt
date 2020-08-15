add_executable(yshapeview  yshapeview.cpp)

set_target_properties(yshapeview  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(yshapeview  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yshapeview  yocto yocto_gui)

