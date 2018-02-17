# Notes on future improvements of Yocto/GL

This file contains notes on future improvements of Yocto.
Please consider this to be just development notes and not any real planning.

## Next

- constexpr
- setup tests
    - regenerate tests one scene at a time
        - report error on scene generation
    - add quiet mode to applications
    - report fatal errors in commands
    - image compression with pngquant
    - lightweight image check with command line diff
    - make results file by filtering applications output
    - output filters applications output
    - consider a tool to run tests

## Deployment

- Doxygen and Sphynx
- insert SHA and version number
- consider amalgamation for yocto
- consider amalgamation for ext

## Obj

- testure info use <</>> operators
    - print before texture path
    - parse up to last argument
- faster OBJ
    - stop using iostreams (darn it)
    - output using fmt or copy its methods?
    - maintain a << API but use custom parsing code

## Tests

- make test suite script
    - consider fast tests and long tests
    - can we test easily yView?

## Scenes

- setup scene repo
- begin importing scenes
- make 4 scene variants
    - original
    - fixed
    - converted OBJ
    - converted GLTF
- consider putting OBJ extensions into its own files?

## Shade uses render buffers

- implement a framebuffer
- hardcode textures inside it
    - create empty texture functions
    - add depth texture for depth
- tutorial at
    http:www.opengl-tutorial.org/intermediate-tutorials/tutorial-14-render-to-texture/

## Tesselation

- update convert functions to new api (?)
- cleanup tesselation in shape
    - remove tesselate once
    - tesselation uses only internal levels

## Math

- frame inverse with flag
- spherical/cartesian conversion
- make stronger the assumption on the use of frames
    - frame inversion
    - documentation
- check random shuffle
- check random number generation for float/double
- check rotation and decompoaition of rotations
   - see euclideanspace.com

## Scene

- share texture info accross GPU/tracer/scene
- make texture info more complete with mirroring and mipmapping
- envmap along z

## Image

- tonemap params to put everywhere
- consider other tone reproduction code
- maybe: make image a simple structure
    - get_pixel, make_image
- remove constructors and accessors from vec/mat/frame

## Low-level code

- serialization with visitor
    - decide if exposing json is reasonable
      for now this is just a matter of compilation time
      later it is best to use a variant type

## Ui

- add angle semantic
- add rotation
- add frame editing with decomposition
- add labels 2,3,4

## Maybe


## Trace

- distributions:
    - move to binary function
    - consider adding an object
    - add a distribution for lights
- cleanup sampling functions everywhere
    - probably removing sample_points/lines/triangles
    - cleanup sampling in ray tracing
    - make lights with single shapes in trace
- add radius in offsetting rays
- simplify trace_point
    - double sided in material functions
    - opacity in material functions
- simplify trace_light
    - maybe include shape directly?
- remove background from point?
- sample background to sum all environments
- envmap sampling
- sobol and cmjs

## BVH

- simplify build functions: can we avoid preallocating nodes?
- maybe put axis with internal
- simplify partition and nth_element function
    - include wrapper functions

## Simple denoiser

- joint bilateral denoiser
- non-local means denoiser

## Apps

- yitrace: check editing
- yitrace: consider update
