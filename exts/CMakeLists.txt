add_subdirectory(fast_float)
add_subdirectory(nlohmann_json)
add_subdirectory(stb_image)
add_subdirectory(tinyexr)
add_subdirectory(cgltf)

if(YOCTO_OPENGL)

set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL " " FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL " " FORCE)
set(GLFW_BUILD_DOCS OFF CACHE BOOL " " FORCE)
set(GLFW_INSTALL OFF CACHE BOOL " " FORCE)
add_subdirectory(glfw)

add_subdirectory(glad)

add_subdirectory(imgui)

endif(YOCTO_OPENGL)

if(YOCTO_EMBREE)
add_subdirectory(embree)
endif(YOCTO_EMBREE)

if(YOCTO_DENOISE)
add_subdirectory(openimagedenoise)
endif(YOCTO_DENOISE)

if(YOCTO_CUDA)
add_subdirectory(embed_ptx)
endif(YOCTO_CUDA)

