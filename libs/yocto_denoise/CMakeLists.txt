if(YOCTO_DENOISE)
  add_definitions(-DYOCTO_DENOISE)
  option(OIDN_APPS "" OFF)
  add_subdirectory("ext/oidn")
  
  add_library(yocto_denoise
    yocto_denoise.cpp
  )

  target_link_libraries(yocto_denoise OpenImageDenoise)
  target_include_directories(yocto_denoise PRIVATE .)
  target_include_directories(yocto_denoise PUBLIC ${CMAKE_SOURCE_DIR}/libs)
  target_compile_definitions(yocto_denoise PRIVATE ext/oidn/include)
  set_target_properties(yocto_denoise PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
endif(YOCTO_DENOISE)

# warning flags
if(APPLE)
  target_compile_options(yocto_gui PRIVATE -Wall)
endif(APPLE)
if(MSVC)
  target_compile_options(yocto_gui PRIVATE /D_CRT_SECURE_NO_WARNINGS /wd4018 /wd4244 /wd4305 /wd4800 /wd4267)
  target_compile_options(yocto_gui PUBLIC /EHsc)
endif(MSVC)