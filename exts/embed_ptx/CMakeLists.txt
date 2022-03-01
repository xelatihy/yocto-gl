project(optix)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} PARENT_SCOPE)

include(embed_ptx)

enable_language(CUDA)

find_package(OptiX)

if (DEFINED ENV{OptiX_INSTALL_DIR})
  find_path(OptiX_ROOT_DIR NAMES include/optix.h PATHS $ENV{OptiX_INSTALL_DIR})
else()
  find_path(OptiX_ROOT_DIR NAMES include/optix.h)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OptiX
  FOUND_VAR OptiX_FOUND
  REQUIRED_VARS
    OptiX_ROOT_DIR
)

# add_library(OptiX::OptiX INTERFACE)
# target_include_directories(OptiX::OptiX INTERFACE ${OptiX_ROOT_DIR}/include)
