cmake_minimum_required (VERSION 3.9)

project (yocto-gl)

option(YOCTO_OPENGL "Build OpenGL apps" ON)
option(YOCTO_EMBREE "Use Intel's Embree raytracer" OFF)
option(YOCTO_TESTING "Enable testing" ON)

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
add_executable(yobjproc apps/yobjproc.cpp)

target_link_libraries(ytrace ygl)
target_link_libraries(yscnproc ygl)
target_link_libraries(yimproc ygl)
target_link_libraries(yobjproc ygl)

if(YOCTO_OPENGL)
    find_package(OpenGL REQUIRED)
    if(APPLE)
        include_directories(/usr/local/include)
        link_directories(/usr/local/lib)
        find_library(GLFW_LIBRARY NAMES glfw3 glfw PATHS /usr/local/lib)
        set(GL_EXTRA_LIBRARIES)
    endif(APPLE)
    if(WIN32)
        include_directories(${CMAKE_SOURCE_DIR}/apps/ext/glfw_w32/include)
        link_directories(${CMAKE_SOURCE_DIR}/apps/ext/glfw_w32/lib-vc2015)
        find_library(GLFW_LIBRARY NAMES glfw3dll PATHS ${CMAKE_SOURCE_DIR}/apps/ext/glfw_w32/lib-vc2015)
        set(GL_EXTRA_LIBRARIES)
    endif(WIN32)
    if(UNIX AND NOT APPLE)
        include_directories(/usr/include /usr/local/include)
        find_library(GLFW_LIBRARY NAMES glfw3 glfw PATHS /usr/lib /usr/local/lib64 /usr/lib64 /usr/local/lib /usr/lib/x86_64-linux-gnu)
        set(GL_EXTRA_LIBRARIES dl)
    endif(UNIX AND NOT APPLE)
    include_directories(apps/)
    add_definitions(-DIMGUI_IMPL_OPENGL_LOADER_GLAD)

    add_library(yopengl apps/yglutils.h apps/yglutils.cpp 
        apps/ext/imgui/imgui.cpp apps/ext/imgui/imgui_draw.cpp apps/ext/imgui/imgui_widgets.cpp
        apps/ext/imgui/imgui_impl_glfw.cpp apps/ext/imgui/imgui_impl_opengl3.cpp apps/ext/glad/glad.c)

    add_executable(yview apps/yview.cpp)
    add_executable(yitrace apps/yitrace.cpp)
    add_executable(yimview apps/yimview.cpp)

    target_link_libraries(yview ygl yopengl ${OPENGL_gl_LIBRARY} ${GLFW_LIBRARY} ${GL_EXTRA_LIBRARIES})
    target_link_libraries(yitrace ygl yopengl ${OPENGL_gl_LIBRARY} ${GLFW_LIBRARY} ${GL_EXTRA_LIBRARIES})
    target_link_libraries(yimview ygl yopengl ${OPENGL_gl_LIBRARY} ${GLFW_LIBRARY} ${GL_EXTRA_LIBRARIES})
endif(YOCTO_OPENGL)

if(YOCTO_TESTING)
enable_testing()
add_executable(yimdiff apps/yimdiff.cpp)
target_link_libraries(yimdiff ygl)
set(test-yimproc ${CMAKE_SOURCE_DIR}/tests/run-tests/test-yimproc.sh)
set(test-ytrace ${CMAKE_SOURCE_DIR}/tests/run-tests/test-ytrace.sh)
set(test-yscnproc ${CMAKE_SOURCE_DIR}/tests/run-tests/test-yscnproc.sh)
set(testdir ${CMAKE_SOURCE_DIR}/tests)
set(outdir ${CMAKE_SOURCE_DIR}/tests/run-tests/output)
set(resdir ${CMAKE_SOURCE_DIR}/tests/run-tests/result)
set(diffdir ${CMAKE_SOURCE_DIR}/tests/run-tests/difference)
foreach(filename "uvgrid.png")
    get_filename_component(name ${filename} NAME_WE)
    add_test(NAME yimproc-${filename}
             COMMAND ${test-yimproc} ${testdir}/textures/${filename} ${outdir}/yimproc-${name}.png ${resdir}/yimproc-${name}.png ${diffdir}/yimproc-${name}.png
             WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})    
endforeach()
foreach(tracer "eyelight" "path" "direct" "path_naive" "path_nomis" "direct_nomis")
    foreach(filename "simple-al.json" "simple-el.json" "cornellbox-al.json" "simplemat-al.json" "simplemat-el.json" "subdiv-al.json"  "subdivfv-al.json" "multimat-al.json" "hairball-al.json" "models-al.json" "emptyelems-al.json")
        get_filename_component(name ${filename} NAME_WE)
        add_test(NAME ytrace-${tracer}-${filename}
                COMMAND ${test-ytrace} ${testdir}/${filename} ${outdir}/ytrace-${tracer}-${name}.png ${resdir}/ytrace-${tracer}-${name}.png ${diffdir}/ytrace-${tracer}-${name}.png "-t ${tracer} -s 4 -r 360"
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    endforeach()
endforeach(tracer)
foreach(format "json" "ybin" "obj")
    foreach(filename "simple-al.json" "simple-el.json" "cornellbox-al.json" "simplemat-al.json" "simplemat-el.json" "subdiv-al.json" "subdivfv-al.json" "multimat-al.json" "models-al.json" "emptyelems-al.json")
        get_filename_component(name ${filename} NAME_WE)
        add_test(NAME yscnproc-${format}-${filename}
                COMMAND ${test-yscnproc} ${testdir}/${filename} ${outdir}/${format}/${name} ${outdir}/${format}/${name}/${name}.${format} ${outdir}/yscnproc-${format}-${name}.png ${resdir}/yscnproc-${format}-${name}.png ${diffdir}/yscnproc-${format}-${name}.png "-t path -s 4 -r 360"
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    endforeach()
endforeach(format)
if(YOCTO_EMBREE)
foreach(tracer "eyelight" "path")
    foreach(filename "simple-al.json" "simple-el.json" "cornellbox-al.json" "simplemat-al.json" "simplemat-el.json" "subdiv-al.json" "subdivfv-al.json" "multimat-al.json" "hairball-al.json" "models-al.json"  "emptyelems-al.json")
        get_filename_component(name ${filename} NAME_WE)
        add_test(NAME ytrace-embree-${tracer}-${filename}
                COMMAND ${test-ytrace} ${testdir}/${filename} ${outdir}/ytrace-embree-${tracer}-${name}.png ${resdir}/ytrace-embree-${tracer}-${name}.png ${diffdir}/ytrace-embree-${tracer}-${name}.png "-t ${tracer} --embree -s 4 -r 360"
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    endforeach()
endforeach(tracer)
endif(YOCTO_EMBREE)
endif(YOCTO_TESTING)
