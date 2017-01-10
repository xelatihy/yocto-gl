# Example applications for YOCTO/GL.

For a complete list of example applications, see the [main readme](../readme.md).

The simplest way to compile the applications is by using `apps/build.osx.sh`,
`apps/build.windows.sh` or `apps/build.linux.sh`. But we also provide a CMake
script as an alternative.

Interactive applications use OpenGL, GLFW and GLEW.
On OsX install the libraries using homebrew.
On Linux install them using you favority package manager
(we tested with Ubuntu/apt-get).
On Windows, we included binary files for both libraries.
