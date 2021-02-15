add_subdirectory(yimage)
add_subdirectory(ytrace)
add_subdirectory(yshape)
add_subdirectory(yscene)

if(YOCTO_OPENGL)
add_subdirectory(yshade)
endif(YOCTO_OPENGL)
