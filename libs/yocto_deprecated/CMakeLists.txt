add_library(yocto_unsupported
  yocto_json.h yocto_json.cpp ext/json.hpp
)

set_target_properties(yocto_unsupported PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

target_include_directories(yocto_unsupported PRIVATE ext/)
target_include_directories(yocto_unsupported PUBLIC ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(yocto_unsupported PUBLIC yocto)
