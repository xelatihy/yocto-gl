# Example applications for YOCTO/GL.

The code depends on GLFW v3 on all platforms and GLEW on Window/Linux.
To build, use the platform build files included in this directory calling them
from the root directory, such as

    ./apps/build.osx.sh

We do not provide build scripts for the OpenGL apps for Windows and Linux since
code dependncy issues are harder to deal with in those systems. You shoul be
able to build them easily though just by looking at the OSX ones.

To make a set of simple tests scenes, use

    ./bin/ytestgen tests

Then you can run `ytrace` for path tracing, `yshade` for quick OpenGL viewing
or `yview` for HDr image viewing. Run the executable with `-h` to get help.
