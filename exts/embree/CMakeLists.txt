add_library(embree SHARED IMPORTED GLOBAL)

if(APPLE)
  if(NOT(CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "arm64"))
    target_include_directories(embree INTERFACE "/usr/local/include")
    set_target_properties(embree PROPERTIES
      IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
      IMPORTED_LOCATION "/usr/local/lib/libembree4.dylib"
      INTERFACE_INCLUDE_DIRECTORIES "/usr/local/include")
  else()
    target_include_directories(embree INTERFACE "/opt/homebrew/include")
    set_target_properties(embree PROPERTIES
      IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
      IMPORTED_LOCATION "/opt/homebrew/lib/libembree4.dylib"
      INTERFACE_INCLUDE_DIRECTORIES "/opt/homebrew/include")
  endif()
endif(APPLE)

if(MSVC)
  target_include_directories(embree INTERFACE "/Program\ Files/Intel/Embree3/include/")
  set_target_properties(embree PROPERTIES
    IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
    IMPORTED_LOCATION "/Program\ Files/Intel/Embree4/bin/embree4.dll"
    IMPORTED_IMPLIB "/Program\ Files/Intel/Embree4/lib/embree4.lib"
    INTERFACE_INCLUDE_DIRECTORIES "/Program\ Files/Intel/Embree4/include/")
endif(MSVC)

if(UNIX AND NOT APPLE)
  target_include_directories(embree INTERFACE /usr/local/include)
  set_target_properties(embree PROPERTIES
    IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
    IMPORTED_LOCATION "/usr/local/lib/libembree4.so"
    INTERFACE_INCLUDE_DIRECTORIES "/usr/local/include")
endif()

add_library(embree::embree ALIAS embree)
