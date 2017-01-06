# Example applications for YOCTO/GL.

1. Make a set of simple tests scenes using `./bin/ytestgen tests`.
2. Render the scene with path tracing with `ytrace` for offline rendering or
   `yitrace` for interactive path tracing.
3. Quickly view your scenes with `yshade`.
4. View rendered images with `yimview` with support for HDR image viewing.
5. Runs rigid body demo with `ysym` for offline simulation or `yisym` for interactive one.

Run the executable with `-h` to get help.

Some of the examples depend on GLFW v3 on all platforms and GLEW on Window/Linux.
The libraries needs to be made available for the build, by changing the paths
in the build files. To build, use the included .sh/.bat files.
