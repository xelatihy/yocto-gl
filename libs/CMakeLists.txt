add_subdirectory(yocto)

if(YOCTO_OPENGL)
add_subdirectory(yocto_gui)
endif(YOCTO_OPENGL)

if(YOCTO_DENOISE)
add_subdirectory(yocto_denoise)
endif(YOCTO_DENOISE)