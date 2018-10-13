# Notes on future improvements of Yocto/GL

This file contains notes on future improvements of Yocto.
Please consider this to be just development notes and not any real planning.

## Next features

- io
    - logging errors but returns always
    - io returns values always and logs errors
    - save functions return booleans
- images
    - need to preserve color space to be able to properly work on images
    - for now, focus on srgb primaries, but start porting other spaces too
    - provide a way to view color spaces in yimview
    - color spaces useful
        - srgb / linear srgb
        - rec2020 / linear rec2020 ?
        - acescg, acescc, acescct
        - some log space
    - practical solution
        - assume srgb, but still label all things clearly
        - do proper srgb coversion
        - imageio supports only srgb for now since the base library do not 
          have access to color profiles
            - check if we can gey access to profiles
        - texture
            - store image4f and image4b data
        - opengl
            - support floating point textures and srgb internal formats 
        - imview
            - this is to be defined
            - problably, do two cases
                - hdr images / linear with srgb primaries
                - ldr image / non-linear with srgb primaries
            - store images in their original form
            - then either provide two separate editing pipelines
            - in general, it feels that we should still have a similar pipeline
              with carefully designed spaces, e.g. using log and always starting
              from a linear representation -> still perform everything in 
              floating point
- yimview
    - add grading
        - check Filmic Worlds implementation
        - check from Unity implementation
        - add tonemapper data
    - add mipmapping
        - easy fix
    - file in
        - tinyDNG
    - add histogram visualization
    - focus on HDR toning, with LDR used as postprocess
    - color correction in both sRGB or after
- color grading
    - schlick formula for bias/gain
    - offset?
    - LUTs
    - tone curve
        - from web?
        - no implementation of Photoshop curve
        - maybe try monotonic spline interpolation
        - otherwise a simple Bezier with constraints
        - use lookup tables for this
- adding quads to shapes
    - integration in sampling
    - integration in intersection
    - integration in interpolation
- cleanup examples
    - bent floor
    - hair
    - crane subdivs
    - utah teapot
    - bunny embed
- BVH
    - add all intersection
    - add cutout to trace
    - simplify build functions
    - maybe put axis with internal
    - simplify partition and nth_element function
        - include wrapper functions

## Cleanup

- generic
    - to_string
    - vec: consider constructors/conversions
    - cmdline with to_string and parse
    - min, lerp, clamp
- path tracing with pbrt MIS
- simple renderers for book
    - simple intersection
    - eyelight
    - anti-aliasing
    - naive path tracing
    - mirrors / glass
    - volume ?
    - product path tracing
- consider introducing trace_point
- refactor renderers into smaller pieces for readability
- direct with deltas
- split apps

## Scene changes

- simpler shapes functions
    - add functions on quads using embree parametrization
    - quad tangent field
    - add quads to shape

## Trace

- clarify light pdfs
- one path
    - compute light pdfs with intersection
- orthonormal basis
    - file:///Users/fabio/Downloads/Duff2017Basis.pdf
- volumes
    - make material models clear
    - start with a simple design
- bump mapping
    - convert bump maps to normal maps internally on load
- displacement
    - apply displacement on tesselation or on triangle vertices directly
- examples
    - hair
    - subdiv
    - displacement
- put back double sided
    - one sided shading
- refraction
    - rough surface refraction
- intersection
    - methods with projections?
        - why shear?
        - triangle
        - line
        - point
    - line/point parametrization that matches pbrt and embree
    - radius in offsetting rays
        - check embree to avoid
- bvh
    - sah
    - opacity
    - embree integration
- see mini path tracer features and include them all
- simple denoiser
    - joint bilateral denoiser
    - non-local means denoiser
    - denoiser from paris

## Port scenes

- bitterli
    - hair
    - render exclude list
    - kitchen
        - bump map
    - car
        - shading bug
        - enable kr
        - metallic paint
    - car2
        - shading bug
        - enable kr
        - metallic paint
    - flip double sided normals on at least large objects
- mcguire
    - exclude list
    - bad textures
    - bad emission
        - lost empire
        - sportscar
    - cameras
    - lights
- pbrt
    - crown material or re-export
    - landscape materials
    - pbrt include parser
- gltf exports

## Additional Notes

### Tonemapping and HDR

- Frostbite HDR
    - from "High Dynamic Range color grading and diksplay in Frostbite"
    - Legacy:
        - hdr/linear -> 
          tonemap/linear -> 
          srgb/non-linear[0-1] -> 
          grading (lut) -> 
          display
    - Current
        - hdr/linear ->
          to lut space/non-linear -> 
          grading (lut) -> 
          to hdr/linear -> 
          tonemap/non-linear for each output (srgb, hdr10) ->
          display
    - Lut space is PQ for now
    - Consider doing it in ACES instead for more compatibility
- COD HDR pipeline
- Filmic Worlds
    - 
- ACES
    - need to read more about this
    - non exact diagram below, but very reasonable
    - raw(non-linear encoding or linear exr) -> [input transform] ->
      aces linear (acescg) -> [look modificartion transform] ->
      aces non-linear (acescc or acescct) -> [color grading] ->
      aces linear (acescg) -> [reference rendering transform] (filminc tone mapping) ->
      aces linear -> [output tranform] ->
      device-dependent color space

### OpenMP

On Apple Clang, you need to add several options to use OpenMP's front end
instead of the standard driver option. This usually looks like
  -Xpreprocessor -fopenmp -lomp

You might need to make sure the lib and include directories are discoverable
if /usr/local is not searched:

  -L/usr/local/opt/libomp/lib -I/usr/local/opt/libomp/include

For CMake, the following flags will cause the OpenMP::OpenMP_CXX target to
be set up correctly:
  -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include" -DOpenMP_CXX_LIB_NAMES="omp" -DOpenMP_omp_LIBRARY=/usr/local/opt/libomp/lib/libomp.dylib
#
## Book

- apps
    - tonemap
    - naivetrace
    - pathtrace
    - itrace
- lib
    - trace.h / trace.cpp
    - traceio.h / traceio.cpp
