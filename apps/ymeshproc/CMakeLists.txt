add_executable(ymeshproc  ymeshproc.cpp)

set_target_properties(ymeshproc  PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(ymeshproc  PRIVATE ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(ymeshproc  yocto)

