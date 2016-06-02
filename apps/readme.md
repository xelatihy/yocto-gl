# Example applications for YOCTO/GL.

The code depends on GLFW v3 on all platforms and GLEW on Window/Linux. These
files are included in the distribution so there should be no need to build
the libraries statically. To build, use either the included 
[Ninja](https://ninja-build.org/) build or the XCode/VisualStudio projects.

To make a set of simple tests scenes, use `./bin/ytestgen tests`.
Then you can run `ytrace` for path tracing, `yshade` for quick OpenGL viewing,
`yview` for HDR image viewing or `ysymrigid` for rigid body simulation. 
Run the executable with `-h` to get help.
