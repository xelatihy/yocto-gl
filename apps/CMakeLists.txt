cmake_minimum_required (VERSION 3.5)

project (yocto-gl)

option(INCLUDE_OPENGL_DEPS "Include OpenGL apps" OFF)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED on)
# set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/../bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# compile flags
if(APPLE)
    add_definitions(-Wno-missing-braces)
endif(APPLE)
if(WIN32)
    add_definitions(/D_CRT_SECURE_NO_WARNINGS /wd4018 /wd4244 /wd4305 /wd4800 /wd4267)
endif(WIN32)

set(app_files
    ../apps/yapp.h
    ../apps/yapp.cpp
    ../apps/tinyply.h
    ../apps/tinyply.cpp
)

if(INCLUDE_OPENGL_DEPS)
    set(app_files
        ${app_files}
        ../apps/yapp_ui.cpp
    )
endif(INCLUDE_OPENGL_DEPS)

add_library(app STATIC ${app_files})

target_link_libraries(app yocto)

add_executable(ysym ../apps/ysym.cpp)
add_executable(ytestgen ../apps/ytestgen.cpp)
add_executable(ytrace ../apps/ytrace.cpp)
add_executable(yobj2gltf ../apps/yobj2gltf.cpp)
add_executable(yimproc ../apps/yimproc.cpp)

target_link_libraries(ysym yocto app)
target_link_libraries(ytestgen yocto app)
target_link_libraries(ytrace yocto app)
target_link_libraries(yobj2gltf yocto)
target_link_libraries(yimproc yocto)

if(UNIX AND NOT APPLE)
    set_target_properties(ytrace PROPERTIES COMPILE_FLAGS -pthread LINK_FLAGS -pthread)
endif(UNIX AND NOT APPLE)

if(INCLUDE_OPENGL_DEPS)
    add_executable(yisym ../apps/yisym.cpp)
    add_executable(yshade ../apps/yshade.cpp)
    add_executable(yitrace ../apps/yitrace.cpp)
    add_executable(yimview ../apps/yimview.cpp)

    target_link_libraries(yisym yocto app)
    target_link_libraries(yshade yocto app)
    target_link_libraries(yimview yocto app)
    target_link_libraries(yitrace yocto app)

    if(UNIX AND NOT APPLE)
        set_target_properties(yitrace PROPERTIES COMPILE_FLAGS -pthread LINK_FLAGS -pthread)
    endif(UNIX AND NOT APPLE)
endif(INCLUDE_OPENGL_DEPS)
