add_subdirectory(yocto)
add_subdirectory(yocto_unsupported)

if(YOCTO_OPENGL)
add_subdirectory(yocto_gui)
endif(YOCTO_OPENGL)
