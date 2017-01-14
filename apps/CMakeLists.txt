cmake_minimum_required (VERSION 3.5)

project (yocto)

option(BUILD_OPENGL_APPS "Build OpenGL apps" ON)

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

if(BUILD_OPENGL_APPS)
    # opengl
    find_package(OpenGL REQUIRED)

    # opengl support libaries
    if(APPLE)
        include_directories(/usr/local/include)
        link_directories(/usr/local/lib)
        find_library(GLFW_LIBRARY NAMES glfw3 glfw PATHS /usr/local/lib)
    endif(APPLE)
    if(WIN32)
    	include_directories(${CMAKE_SOURCE_DIR}/w32/include)
        link_directories(${CMAKE_SOURCE_DIR}/w32/lib/lib-vc2015)
        find_library(GLEW_LIBRARIES NAMES glew32 PATHS ${CMAKE_SOURCE_DIR}/w32/lib-vc2015)
        find_library(GLFW_LIBRARY NAMES glfw3dll PATHS ${CMAKE_SOURCE_DIR}/w32/lib-vc2015)
    endif(WIN32)
    if(UNIX AND NOT APPLE)
        include_directories(/usr/include /usr/local/include)
        find_library(GLFW_LIBRARY NAMES glfw3 glfw PATHS /usr/lib /usr/local/lib64 /usr/lib64 /usr/local/lib /usr/lib/x86_64-linux-gnu)
        find_package(GLEW REQUIRED)
    endif(UNIX AND NOT APPLE)
endif(BUILD_OPENGL_APPS)

set(yocto_files
    ../yocto/yocto_cmd.h
    ../yocto/yocto_shape.h
    ../yocto/yocto_obj.h
    ../yocto/yocto_obj.cpp
    ../yocto/yocto_gltf.h
    ../yocto/yocto_gltf.cpp
    ../yocto/yocto_math.h
    ../yocto/yocto_trace.h
    ../yocto/yocto_trace.cpp
    ../yocto/yocto_sym.h
    ../yocto/yocto_sym.cpp
    ../yocto/yocto_bvh.h
    ../yocto/yocto_bvh.cpp
    ../yocto/yocto_shape.h
    ../yocto/yocto_shape.cpp
)

set(app_files
    ../apps/yapp.h
    ../apps/yapp.cpp
)

if(BUILD_OPENGL_APPS)
    set(yocto_files
        ${yocto_files}
        ../yocto/yocto_glu.h
        ../yocto/yocto_glu.cpp
    )
    set(app_files
        ${app_files}
        ../apps/yui.h
        ../apps/yui.cpp
    )
endif(BUILD_OPENGL_APPS)

add_library(yocto STATIC ${yocto_files})
add_library(app STATIC ${app_files})

target_link_libraries(app yocto)

add_executable(ysym ../apps/ysym.cpp)
add_executable(ytestgen ../apps/ytestgen.cpp)
add_executable(ytrace ../apps/ytrace.cpp)
add_executable(yobj2gltf ../apps/yobj2gltf.cpp)

target_link_libraries(ysym yocto app)
target_link_libraries(ytestgen yocto app)
target_link_libraries(ytrace yocto app)
target_link_libraries(yobj2gltf yocto)

if(UNIX AND NOT APPLE)
    set_target_properties(ytrace PROPERTIES COMPILE_FLAGS -pthread LINK_FLAGS -pthread)
endif(UNIX AND NOT APPLE)

if(BUILD_OPENGL_APPS)
    add_executable(yisym ../apps/yisym.cpp ../apps/yapp_ui.cpp ${yall_glfiles})
    add_executable(yshade ../apps/yshade.cpp ../apps/yapp_ui.cpp ${yall_glfiles})
    add_executable(yitrace ../apps/yitrace.cpp ../apps/yapp_ui.cpp ${yall_glfiles})
    add_executable(yimview ../apps/yimview.cpp ../apps/yui.h ../apps/yui.cpp ../yocto/yocto_cmd.h ../yocto/yocto_math.h ../yocto/yocto_glu.h ../yocto/yocto_glu.cpp)

    target_link_libraries(yisym yocto app ${OPENGL_gl_LIBRARY} ${GLFW_LIBRARY} ${GLEW_LIBRARIES})
    target_link_libraries(yshade yocto app ${OPENGL_gl_LIBRARY} ${GLFW_LIBRARY} ${GLEW_LIBRARIES})
    target_link_libraries(yimview yocto app ${OPENGL_gl_LIBRARY} ${GLFW_LIBRARY} ${GLEW_LIBRARIES})
    target_link_libraries(yitrace yocto app ${OPENGL_gl_LIBRARY} ${GLFW_LIBRARY} ${GLEW_LIBRARIES})

    if(UNIX AND NOT APPLE)
        set_target_properties(yitrace PROPERTIES COMPILE_FLAGS -pthread LINK_FLAGS -pthread)
    endif(UNIX AND NOT APPLE)
endif(BUILD_OPENGL_APPS)
