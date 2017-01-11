# Example applications for YOCTO/GL.

For a complete list of example applications, see the [main readme](../readme.md).

The simplest way to compile the non-interative applications is by using
`apps/build.osx.sh`, `apps/build.windows.sh` or `apps/build.linux.sh`. These
applications have no external dependency beside the C++ STL.

For both the interactive and non-interactive applicaitons, we provide a
CMake script. Interactive applications use OpenGL, GLFW and GLEW.
On OSX install the libraries using homebrew. On Linux install them using you
favorite package manager (we tested with Ubuntu/apt-get).
On Windows, we included binary files for both libraries.
