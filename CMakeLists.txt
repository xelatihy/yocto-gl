cmake_minimum_required (VERSION 3.5)

project (yocto-gl)

option(YOCTO_OPENGL "Build OpenGL apps" ON)
option(YOCTO_EMBREE "Use Intel's Embree raytracer" ON)

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
if(YOCTO_OPENGL)
    add_definitions(-DYGL_OPENGL=1)
else(YOCTO_OPENGL)
    add_definitions(-DYGL_OPENGL=0)
endif(YOCTO_OPENGL)
if(YOCTO_EMBREE)
    add_definitions(-DYGL_EMBREE=1)
else(YOCTO_EMBREE)
    add_definitions(-DYGL_EMBREE=0)
endif(YOCTO_EMBREE)

# platform compile flags
if(APPLE)
    add_definitions(-Wall)
    if(YOCTO_EMBREE)
    include_directories(/opt/local/include)
    link_directories(/opt/local/lib)
    endif(YOCTO_EMBREE)
endif(APPLE)
if(MSVC)
    add_definitions(/D_CRT_SECURE_NO_WARNINGS /wd4018 /wd4244 /wd4305 /wd4800 /wd4267)
    set(CMAKE_CXX_FLAGS "/EHsc")
endif(MSVC)

add_library(ygl 
    yocto/ygl.h yocto/ygl.cpp
    yocto/yglio.h yocto/yglio.cpp
    yocto/yobj.h yocto/ygltf.h
)

if(UNIX AND NOT APPLE)
    find_package(Threads REQUIRED)
    target_link_libraries(ygl Threads::Threads)
endif(UNIX AND NOT APPLE)

if(YOCTO_EMBREE)
target_link_libraries(ygl libembree3.dylib)
endif(YOCTO_EMBREE)

add_executable(ytrace apps/ytrace.cpp)
add_executable(yscnproc apps/yscnproc.cpp)
add_executable(yimproc apps/yimproc.cpp)
add_executable(ytonemap apps/ytonemap.cpp)
add_executable(ytrace2 apps/ytrace2.cpp)
add_executable(yobjproc apps/yobjproc.cpp)

target_link_libraries(ytrace ygl)
target_link_libraries(yscnproc ygl)
target_link_libraries(yimproc ygl)
target_link_libraries(ytonemap ygl)
target_link_libraries(ytrace2 ygl)
target_link_libraries(yobjproc ygl)

if(YOCTO_OPENGL)
    find_package(OpenGL REQUIRED)
    if(APPLE)
        include_directories(/usr/local/include)
        link_directories(/usr/local/lib)
        find_library(GLFW_LIBRARY NAMES glfw3 glfw PATHS /usr/local/lib)
        add_definitions(-DIMGUI_IMPL_OPENGL_LOADER_CUSTOM=<OpenGL/gl3.h>)
    endif(APPLE)
    if(WIN32)
        include_directories(${CMAKE_SOURCE_DIR}/apps/glw32/include)
        link_directories(${CMAKE_SOURCE_DIR}/apps/glw32/lib-vc2015)
        find_library(GLEW_LIBRARIES NAMES glew32 PATHS ${CMAKE_SOURCE_DIR}/apps/glw32/lib-vc2015)
        find_library(GLFW_LIBRARY NAMES glfw3dll PATHS ${CMAKE_SOURCE_DIR}/apps/glw32/lib-vc2015)
    endif(WIN32)
    if(UNIX AND NOT APPLE)
        include_directories(/usr/include /usr/local/include)
        find_library(GLFW_LIBRARY NAMES glfw3 glfw PATHS /usr/lib /usr/local/lib64 /usr/lib64 /usr/local/lib /usr/lib/x86_64-linux-gnu)
        find_package(GLEW REQUIRED)
    endif(UNIX AND NOT APPLE)

    add_executable(yview apps/yview.cpp       
        apps/imgui/imgui.cpp apps/imgui/imgui_draw.cpp apps/imgui/imgui_widgets.cpp
        apps/imgui/imgui_impl_glfw.cpp apps/imgui/imgui_impl_opengl3.cpp)
    add_executable(yitrace apps/yitrace.cpp        
        apps/imgui/imgui.cpp apps/imgui/imgui_draw.cpp apps/imgui/imgui_widgets.cpp
        apps/imgui/imgui_impl_glfw.cpp apps/imgui/imgui_impl_opengl2.cpp)
    add_executable(yimview apps/yimview.cpp        
        apps/imgui/imgui.cpp apps/imgui/imgui_draw.cpp apps/imgui/imgui_widgets.cpp
        apps/imgui/imgui_impl_glfw.cpp apps/imgui/imgui_impl_opengl2.cpp)
    add_executable(yimview2 apps/yimview2.cpp)
    add_executable(yitrace2 apps/yitrace2.cpp)

    target_link_libraries(yview ygl ${OPENGL_gl_LIBRARY} ${GLFW_LIBRARY} ${GLEW_LIBRARIES})
    target_link_libraries(yitrace ygl ${OPENGL_gl_LIBRARY} ${GLFW_LIBRARY} ${GLEW_LIBRARIES})
    target_link_libraries(yimview ygl ${OPENGL_gl_LIBRARY} ${GLFW_LIBRARY} ${GLEW_LIBRARIES})
    target_link_libraries(yimview2 ygl ${OPENGL_gl_LIBRARY} ${GLFW_LIBRARY})
    target_link_libraries(yitrace2 ygl ${OPENGL_gl_LIBRARY} ${GLFW_LIBRARY})
endif(YOCTO_OPENGL)
