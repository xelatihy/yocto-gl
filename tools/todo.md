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

- api for trace with callbacks
- review sample_light_pdf -> explicit normal to be carried around
    - or remove sample_light and use sample_direction
- intersection in scene with other name
    - maybe carry environment here?
- sampling shapes
    - remove sample_points, sample_lines, etc...
    - support uniform sampling in sample_shape
- deltas with only one function call instead of split
- split eval functions from point
    - eval emission
        - simplify?
    - eval environment
        - flip direction?
- redirect log to imgui
- BUG: double sided for kt
- BUG: opacity
- BUG: add default material on load
- opacity in intersection
- sample pdfs based on intersection record
    - probability based on intersection point (pos, norm, ei)
        - fold ei into point
- environment based on only w?
    - do not need point for the environment any more?
    - Le -> emission
    - ei -> pdf
- trace point
    - should we have eval_pos, eval_norm, eval_emission, eval_brdf?
        - only problem is where to put normal map
    - use intersection when needed
        - seems to possible have better semantic for envmap?
    - remvoe point out of the equation
- intersection return projections for lines and points
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
- add radius in offsetting rays
- sample background to sum all environments
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
