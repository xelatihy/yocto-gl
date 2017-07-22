cmake_minimum_required (VERSION 3.5)

project (yocto-gl)

option(YOCTO_NO_OPENGL "Skip OpenGL apps" OFF)
option(YOCTO_CPP11 "Compile as C++11" OFF)

if(NOT YOCTO_CPP11)
set(CMAKE_CXX_STANDARD 14)
else(NOT YOCTO_CPP11)
set(CMAKE_CXX_STANDARD 11)
endif(NOT YOCTO_CPP11)

set(CMAKE_CXX_STANDARD_REQUIRED on)
# set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
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

if(NOT YOCTO_NO_OPENGL)
    # opengl
    find_package(OpenGL REQUIRED)

    # opengl support libaries
    if(APPLE)
        include_directories(/usr/local/include)
        link_directories(/usr/local/lib)
        find_library(GLFW_LIBRARY NAMES glfw3 glfw PATHS /usr/local/lib)
    endif(APPLE)
    if(WIN32)
    	include_directories(${CMAKE_SOURCE_DIR}/apps/w32/include)
        link_directories(${CMAKE_SOURCE_DIR}/apps/w32/lib/lib-vc2015)
        find_library(GLEW_LIBRARIES NAMES glew32 PATHS ${CMAKE_SOURCE_DIR}/apps/w32/lib-vc2015)
        find_library(GLFW_LIBRARY NAMES glfw3dll PATHS ${CMAKE_SOURCE_DIR}/apps/w32/lib-vc2015)
    endif(WIN32)
    if(UNIX AND NOT APPLE)
        include_directories(/usr/include /usr/local/include)
        find_library(GLFW_LIBRARY NAMES glfw3 glfw PATHS /usr/lib /usr/local/lib64 /usr/lib64 /usr/local/lib /usr/lib/x86_64-linux-gnu)
        find_package(GLEW REQUIRED)
    endif(UNIX AND NOT APPLE)
endif(NOT YOCTO_NO_OPENGL)

set(yocto_files
    yocto_math.h
    yocto_utils.h
    yocto_obj.h
    yocto_obj.cpp
    yocto_gltf.h
    yocto_gltf.cpp
    yocto_trace.h
    yocto_trace.cpp
    yocto_sym.h
    yocto_sym.cpp
    yocto_bvh.h
    yocto_bvh.cpp
    yocto_img.h
    yocto_img.cpp
)

if(NOT YOCTO_NO_OPENGL)
    set(yocto_files
        ${yocto_files}
        yocto_glu.h
        yocto_glu.cpp
        yocto_gui.h
        yocto_gui.cpp
        ext/imgui/imgui.cpp
        ext/imgui/imgui_draw.cpp
        ext/imgui/imgui_demo.cpp
        ext/imgui/imgui_impl_glfw_gl3.cpp
    )
    if(NOT APPLE)
    set(yocto_files
        ${yocto_files}
        ext/glad/glad.h
        ext/glad/khrplatform.h
        ext/glad/glad.c
    )
    endif(NOT APPLE)

endif(NOT YOCTO_NO_OPENGL)

add_library(yocto STATIC ${yocto_files})

if(UNIX AND NOT APPLE)
    find_package(Threads REQUIRED)
    target_link_libraries(yocto Threads::Threads)
endif(UNIX AND NOT APPLE)

if(NOT YOCTO_NO_OPENGL)
    target_link_libraries(yocto ${OPENGL_gl_LIBRARY} ${GLFW_LIBRARY} ${GLEW_LIBRARIES})
endif(NOT YOCTO_NO_OPENGL)
