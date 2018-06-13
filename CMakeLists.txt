cmake_minimum_required (VERSION 3.5)

project (yocto-gl)

option(YOCTO_OPENGL "Build OpenGL apps" ON)
option(YOCTO_EXPERIMENTAL "Build experimental apps" OFF)
option(YOCTO_TOOLS "Build tools apps" OFF)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# compile flags
if(APPLE)
    add_definitions(-Wall)
endif(APPLE)
if(MSVC)
    add_definitions(/D_CRT_SECURE_NO_WARNINGS /wd4018 /wd4244 /wd4305 /wd4800 /wd4267)
    set(CMAKE_CXX_FLAGS "/EHsc")
endif(MSVC)
if(YOCTO_OPENGL)
    add_definitions(-DYGL_OPENGL=1)
else(YOCTO_OPENGL)
    add_definitions(-DYGL_OPENGL=0)
endif(YOCTO_OPENGL)

set(YOCTO_FILES
    yocto/yocto_gl.h yocto/yocto_gl.cpp
    yocto/yocto_glio.h yocto/yocto_glio.cpp
    yocto/yocto_obj.h yocto/yocto_gltf.h
)

if(YOCTO_OPENGL)
    find_package(OpenGL REQUIRED)
    if(APPLE)
        include_directories(/usr/local/include)
        link_directories(/usr/local/lib)
        find_library(GLFW_LIBRARY NAMES glfw3 glfw PATHS /usr/local/lib)
    endif(APPLE)
    if(WIN32)
        include_directories(${CMAKE_SOURCE_DIR}/apps/w32/include)
        link_directories(${CMAKE_SOURCE_DIR}/apps/w32/lib-vc2015)
        find_library(GLEW_LIBRARIES NAMES glew32 PATHS ${CMAKE_SOURCE_DIR}/apps/w32/lib-vc2015)
        find_library(GLFW_LIBRARY NAMES glfw3dll PATHS ${CMAKE_SOURCE_DIR}/apps/w32/lib-vc2015)
    endif(WIN32)
    if(UNIX AND NOT APPLE)
        include_directories(/usr/include /usr/local/include)
        find_library(GLFW_LIBRARY NAMES glfw3 glfw PATHS /usr/lib /usr/local/lib64 /usr/lib64 /usr/local/lib /usr/lib/x86_64-linux-gnu)
        find_package(GLEW REQUIRED)
    endif(UNIX AND NOT APPLE)
    set(YOCTO_FILES 
        ${YOCTO_FILES}
        yocto/yocto_glutils.h yocto/yocto_glutils.cpp        
        yocto/ext/imgui/imgui.cpp yocto/ext/imgui/imgui_draw.cpp yocto/ext/imgui/imgui_demo.cpp yocto/ext/imgui/imgui_impl_glfw_gl3.cpp yocto/ext/imgui/imgui_extra_fonts.cpp    )
endif(YOCTO_OPENGL)

add_library(yocto_gl ${YOCTO_FILES})

if(YOCTO_OPENGL)
target_link_libraries(yocto_gl ${OPENGL_gl_LIBRARY} ${GLFW_LIBRARY} ${GLEW_LIBRARIES})
endif(YOCTO_OPENGL)

if(UNIX AND NOT APPLE)
    find_package(Threads REQUIRED)
    target_link_libraries(yocto_gl Threads::Threads)
endif(UNIX AND NOT APPLE)

add_executable(ytrace apps/ytrace.cpp)
add_executable(yscnproc apps/yscnproc.cpp)
add_executable(yimproc apps/yimproc.cpp)

target_link_libraries(ytrace yocto_gl)
target_link_libraries(yscnproc yocto_gl)
target_link_libraries(yimproc yocto_gl)

if(YOCTO_OPENGL)
    add_executable(yview apps/yview.cpp apps/yapp_ui.h  apps/yapp_ui.cpp)
    add_executable(yitrace apps/yitrace.cpp apps/yapp_ui.h  apps/yapp_ui.cpp)
    add_executable(yimview apps/yimview.cpp)

    target_link_libraries(yview yocto_gl)
    target_link_libraries(yitrace yocto_gl)
    target_link_libraries(yimview yocto_gl)
endif(YOCTO_OPENGL)

if(YOCTO_TOOLS)
    # add_executable(ytungsten tools/converters/ytungsten.cpp)
    add_executable(ypbrt tools/converters/ypbrt.cpp)

    # target_link_libraries(ytungsten yocto_gl)
    target_link_libraries(ypbrt yocto_gl)
endif(YOCTO_TOOLS)
