cmake_minimum_required (VERSION 3.9)

project (yocto-gl)

option(YOCTO_OPENGL "Build OpenGL apps" ON)
option(YOCTO_EMBREE "Use Intel's Embree raytracer" OFF)
option(YOCTO_TESTING "Enable testing" ON)

set(CMAKE_CXX_STANDARD 17)
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
    add_definitions(-DYOCTO_OPENGL=1)
else(YOCTO_OPENGL)
    add_definitions(-DYOCTO_OPENGL=0)
endif(YOCTO_OPENGL)
if(YOCTO_EMBREE)
    add_definitions(-DYOCTO_EMBREE=1)
else(YOCTO_EMBREE)
    add_definitions(-DYOCTO_EMBREE=0)
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

add_library(yocto 
    yocto/yocto_math.h yocto/yocto_random.h
    yocto/yocto_bvh.h yocto/yocto_bvh.cpp
    yocto/yocto_shape.h yocto/yocto_shape.cpp 
    yocto/yocto_image.h yocto/yocto_image.cpp
    yocto/yocto_scene.h yocto/yocto_scene.cpp
    yocto/yocto_trace.h yocto/yocto_trace.cpp
    yocto/yocto_obj.h yocto/yocto_obj.cpp 
    yocto/yocto_pbrt.h yocto/yocto_pbrt.cpp 
    yocto/yocto_sceneio.h yocto/yocto_sceneio.cpp
)

if(UNIX AND NOT APPLE)
    find_package(Threads REQUIRED)
    target_link_libraries(yocto Threads::Threads)
endif(UNIX AND NOT APPLE)

if(YOCTO_EMBREE)
target_link_libraries(yocto libembree3.dylib)
endif(YOCTO_EMBREE)

add_executable(yscntrace apps/yscntrace.cpp)
add_executable(yscnproc apps/yscnproc.cpp)
add_executable(yscnmerge apps/yscnmerge.cpp)
add_executable(yimgproc apps/yimgproc.cpp)
add_executable(ymshproc apps/ymshproc.cpp)
add_executable(yislandproc apps/yislandproc.cpp)

target_link_libraries(yscntrace yocto)
target_link_libraries(yscnproc yocto)
target_link_libraries(yscnmerge yocto)
target_link_libraries(yimgproc yocto)
target_link_libraries(ymshproc yocto)
target_link_libraries(yislandproc yocto)

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

    add_library(yocto_opengl apps/yocto_opengl.h apps/yocto_opengl.cpp 
        apps/ext/imgui/imgui.cpp apps/ext/imgui/imgui_draw.cpp apps/ext/imgui/imgui_widgets.cpp
        apps/ext/imgui/imgui_impl_glfw.cpp apps/ext/imgui/imgui_impl_opengl3.cpp apps/ext/glad/glad.c)

    add_executable(yscnview apps/yscnview.cpp)
    add_executable(yscnitrace apps/yscnitrace.cpp)
    add_executable(yimgview apps/yimgview.cpp)

    target_link_libraries(yscnview yocto yocto_opengl ${OPENGL_gl_LIBRARY} ${GLFW_LIBRARY} ${GL_EXTRA_LIBRARIES})
    target_link_libraries(yscnitrace yocto yocto_opengl ${OPENGL_gl_LIBRARY} ${GLFW_LIBRARY} ${GL_EXTRA_LIBRARIES})
    target_link_libraries(yimgview yocto yocto_opengl ${OPENGL_gl_LIBRARY} ${GLFW_LIBRARY} ${GL_EXTRA_LIBRARIES})
endif(YOCTO_OPENGL)

if(YOCTO_TESTING)
enable_testing()
add_executable(yimgdiff apps/yimgdiff.cpp)
target_link_libraries(yimgdiff yocto)
set(test-yimgproc ${CMAKE_SOURCE_DIR}/tests/_results/test-yimgproc.sh)
set(test-yscntrace ${CMAKE_SOURCE_DIR}/tests/_results/test-yscntrace.sh)
set(test-yscnproc ${CMAKE_SOURCE_DIR}/tests/_results/test-yscnproc.sh)
set(testdir ${CMAKE_SOURCE_DIR}/tests)
set(outdir ${CMAKE_SOURCE_DIR}/tests/_output)
set(resdir ${CMAKE_SOURCE_DIR}/tests/_results)
set(diffdir ${CMAKE_SOURCE_DIR}/tests/_difference)
# image io
foreach(filename "uvgrid.png")
    get_filename_component(name ${filename} NAME_WE)
    add_test(NAME yimgproc-${filename}
             COMMAND ${test-yimgproc} ${testdir}/textures/${filename} ${outdir}/yimgproc-${name}.png ${resdir}/yimgproc-${name}.png ${diffdir}/yimgproc-${name}.png
             WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})    
endforeach()
# procedurla images
foreach(name "images1" "images2")
    add_test(NAME yimgproc-${name}
             COMMAND ${test-yimgproc} ::yocto::${name}::/${testdir}/${name}.png ${outdir}/yimgproc-${name}.png ${resdir}/yimgproc-${name}.png ${diffdir}/yimgproc-${name}.png
             WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})    
endforeach()
# eyelight
foreach(filename "features1.yaml" "materials1.yaml" "materials2.yaml" "materials3.yaml" "materials4.yaml" "shapes1.yaml" "shapes2.yaml" "shapes3.yaml" "arealights1.yaml" "environments1.yaml")
    set(tracer "eyelight")
    get_filename_component(name ${filename} NAME_WE)
    add_test(NAME yscntrace-${tracer}-${filename}
            COMMAND ${test-yscntrace} ${testdir}/${filename} ${outdir}/yscntrace-${tracer}-${name}.png ${resdir}/yscntrace-${tracer}-${name}.png ${diffdir}/yscntrace-${tracer}-${name}.png "-t ${tracer} -s 4 -r 640"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endforeach()
# path
foreach(filename "features1.yaml" "materials1.yaml" "materials2.yaml" "materials3.yaml" "materials4.yaml" "shapes1.yaml" "shapes2.yaml" "shapes3.yaml" "arealights1.yaml" "environments1.yaml")
    set(tracer "path")
    get_filename_component(name ${filename} NAME_WE)
    add_test(NAME yscntrace-${tracer}-${filename}
            COMMAND ${test-yscntrace} ${testdir}/${filename} ${outdir}/yscntrace-${tracer}-${name}.png ${resdir}/yscntrace-${tracer}-${name}.png ${diffdir}/yscntrace-${tracer}-${name}.png "-t ${tracer} -s 16 -r 640"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endforeach()
# naive
foreach(filename "features1.yaml" "materials1.yaml" "materials2.yaml" "materials3.yaml" "shapes1.yaml" "shapes2.yaml" "shapes3.yaml" "arealights1.yaml" "environments1.yaml")
    set(tracer "naive")
    get_filename_component(name ${filename} NAME_WE)
    add_test(NAME yscntrace-${tracer}-${filename}
            COMMAND ${test-yscntrace} ${testdir}/${filename} ${outdir}/yscntrace-${tracer}-${name}.png ${resdir}/yscntrace-${tracer}-${name}.png ${diffdir}/yscntrace-${tracer}-${name}.png "-t ${tracer} -s 16 -r 640"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endforeach()
if(YOCTO_EMBREE)
# embree
foreach(filename "features1.yaml" "materials1.yaml" "materials2.yaml" "materials3.yaml" "shapes1.yaml" "shapes2.yaml" "shapes3.yaml" "arealights1.yaml" "environments1.yaml")
    set(tracer "path")
    get_filename_component(name ${filename} NAME_WE)
    add_test(NAME yscntrace-embree-${tracer}-${filename}
            COMMAND ${test-yscntrace} ${testdir}/${filename} ${outdir}/yscntrace-embree-${tracer}-${name}.png ${resdir}/yscntrace-embree-${tracer}-${name}.png ${diffdir}/yscntrace-embree-${tracer}-${name}.png "-t ${tracer} --bvh-embree -s 16 -r 640"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endforeach()
endif(YOCTO_EMBREE)
# yaml output
foreach(filename "features1.yaml" "shapes1.yaml" "shapes2.yaml" "shapes3.yaml")
    set(format "yaml")
    get_filename_component(name ${filename} NAME_WE)
    add_test(NAME yscnproc-${format}-${filename}
            COMMAND ${test-yscnproc} ${testdir}/${filename} ${outdir}/${format}/${name} ${outdir}/${format}/${name}/${name}.${format} ${outdir}/yscnproc-${format}-${name}.png ${resdir}/yscnproc-${format}-${name}.png ${diffdir}/yscnproc-${format}-${name}.png "-t path -s 8 -r 640"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endforeach()
# obj output
foreach(filename "features1.yaml" "shapes1.yaml" "shapes2.yaml" "shapes3.yaml")
    set(format "obj")
    get_filename_component(name ${filename} NAME_WE)
    add_test(NAME yscnproc-${format}-${filename}
            COMMAND ${test-yscnproc} ${testdir}/${filename} ${outdir}/${format}/${name} ${outdir}/${format}/${name}/${name}.${format} ${outdir}/yscnproc-${format}-${name}.png ${resdir}/yscnproc-${format}-${name}.png ${diffdir}/yscnproc-${format}-${name}.png "-t path -s 8 -r 640"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endforeach()
# pbrt output
foreach(filename "features1.yaml" "shapes1.yaml" "shapes2.yaml" "shapes3.yaml")
    set(format "pbrt")
    get_filename_component(name ${filename} NAME_WE)
    add_test(NAME yscnproc-${format}-${filename}
            COMMAND ${test-yscnproc} ${testdir}/${filename} ${outdir}/${format}/${name} ${outdir}/${format}/${name}/${name}.${format} ${outdir}/yscnproc-${format}-${name}.png ${resdir}/yscnproc-${format}-${name}.png ${diffdir}/yscnproc-${format}-${name}.png "-t path -s 8 -r 640"
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endforeach()
endif(YOCTO_TESTING)
