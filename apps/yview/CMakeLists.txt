add_executable(yview  yview.cpp)

set_target_properties(yview  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(yview  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yview PRIVATE yocto)
