# Notes on future improvements of Yocto/GL

This file contains notes on future improvements of Yocto.
Please consider this to be just development notes and not any real planning.

## Main features

- better tests
- better material rendering
- better rendering
- interactive procedural shapes
- prepare for research on procedural components

## Bugs

- lines test
- command line flags for booleans

## Cleanup

- delete
    - compositing if not used

## Tests

- cornell box with light types
- simple scene with pointers
- simple scene with area/env lights

## Trace

- review loops for intersect missed test
    - consider removing environment point
- opacity as lobe
- deltas mixed back in
- matcap rendering
- area light test scene
- bump/normal mapping
- shape trimming in intersection
- double-sided brdf
- fresnel
    - diffuse formula
    - scale with roughness
    - transmission formula
    - fresnel in coefficients
    - fresnel in weights
- highlights are too soft in bitterli scenes
- fast distribution sampling
- add radius in offsetting rays
- sample background to sum all environments
- envmap sampling
    - simplest case, pick pixels as 1D distribution
- cleanup sampling functions everywhere
    - probably removing sample_points/lines/triangles
    - cleanup sampling in ray tracing
    - make lights with single shapes in trace
- simple denoiser
    - joint bilateral denoiser
    - non-local means denoiser
    - denoiser from paris
- evaluate a big simplifying pass

## Tesselated shapes

- add tesselated shape
    - multiple materials
    - face varying here
    - gets a pointer to shape array
    - node can point to this 

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

## Test scenes: simnplify shape generation

- remove some tesselation function
- use subdivs for shader ball
- substance-like shader ball
- bent floor
- 0 roughness

## OpenGL/Trace

- investigate bump map on GPU
    - https://www.opengl.org/discussion_boards/showthread.php/162857-Computing-the-tangent-space-in-the-fragment-shader
    - http://jbit.net/~sparky/sfgrad_bump/mm_sfgrad_bump.pdf

## yView

- double sided in OBJ
- cutout not working on yview
    - consider cutout by default

## BVH

- SAH based build?
- simplify build code: can we avoid preallocating nodes?
- move away from special functions in BVH?
    - always use sort
    - provide a sort buffer
- add cutout to trace
- simplify build functions: 
- maybe put axis with internal
- simplify partition and nth_element function
    - include wrapper functions

## Math

- consider removing const refs
- frame inverse with flag
- check random shuffle
- check random number generation for float/double
