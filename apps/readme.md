# Example applications for YOCTO/GL.

The code depends on GLFW v3 on all platforms and GLEW on Window/Linux. The
libraries needs to be made available for the build, by changing the paths in the build files. To build, use the included .sh/.bat files.

To make a set of simple tests scenes, use `./bin/ytestgen tests`.
Then you can run `ytrace` for path tracing, `yshade` for quick OpenGL viewing,
`yimview` for HDR image viewing or `ysym` for rigid body simulation. 
Run the executable with `-h` to get help.
