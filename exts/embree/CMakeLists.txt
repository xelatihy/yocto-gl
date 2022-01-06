add_library(embree SHARED IMPORTED GLOBAL)

if(APPLE)
  target_include_directories(embree INTERFACE /usr/local/include)
  set_target_properties(embree PROPERTY
    IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
    IMPORTED_LOCATION "/usr/local/lib/libembree3.dylib"
    INTERFACE_INCLUDE_DIRECTORIES "/usr/local/include")
endif(APPLE)
if(MSVC)
  target_include_directories(embree INTERFACE "/Program\ Files/Intel/Embree3/include/")
  set_target_properties(embree PROPERTIES
    IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
    IMPORTED_LOCATION "/Program\ Files/Intel/Embree3/bin/embree3.dll"
    IMPORTED_IMPLIB   "/Program\ Files/Intel/Embree3/lib/embree3.lib"
    INTERFACE_INCLUDE_DIRECTORIES "/Program\ Files/Intel/Embree3/include/")
endif(MSVC)
if(UNIX AND NOT APPLE)
  target_include_directories(embree INTERFACE /usr/local/include)
  set_target_properties(embree PROPERTY
    IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
    IMPORTED_LOCATION "/usr/local/lib/libembree3.so"
    INTERFACE_INCLUDE_DIRECTORIES "/usr/local/include")
endif()

add_library(embree::embree ALIAS embree)
