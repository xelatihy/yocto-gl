# Notes on future improvements of Yocto/GL

This file contains notes on future improvements of Yocto.
Please consider this to be just development notes and not any real planning.

## Main features

- better tests
- better material rendering
- better rendering
- interactive procedural shapes
- prepare for research on procedural components

## Simplifications

- bvh
    - do not sort elements
    - simpler instance api
- bugs
    - eval_pos(env) should apply frame
    - implement eval_pos(ist)
    - trace should use these
- add tesselated shape
    - multiple materials
    - face varying here
    - gets a pointer to shape array
    - node can point to this 

## Trace

- simplification
    - trace functions returns vec4f 
    - double-sided brdf
    - remove none point?
    - is point needed?
        - introduce brdf, pos, norm
    - trace_lights
        - remove
        - add node point
        - solve the stupid problem
        - lights as pointers
        - maybe integrate lights back into the scene
        - add pointer to node if helpful
    - trace_pixel and trace_sample
        - do we need trace_pixels?
    - remove environment map point?
    - sample all lights only
    - back to using pdfs
    - remove distributions
- variants
    - steve's
    - one sample mis only
- lower energy on tests
- test on car paint
- fresnel
    - diffuse formula
    - scale with roughness
    - transmission formula
    - fresnel in coefficients
    - fresnel in weights
- brdf
    - deltas without delta flag
- highlights are too soft in bitterli scenes
- roussian roulette on weight
- samplers
    - sobol
    - cmjs
    - adaptive sampling ala tungsgen
- light sampling
    - possible bug in light weight
    - envmap sampling
    - path trace with explicit light sampling
    - check tungsten light smapling
    - eval_direct function
    - mis in params and not renderer?
- bump mapping and normal mapping
    - bump mapping frame
- add cutout to trace
    - simple scheme: recurse intersect_scene on opacity
    - better scheme: filter intersections
- add transmission to trace
    - better scheme: filter intersection
- distributions
    - is using simple distribution, go back to CDF only
        - sample_distribution
    - move to binary function
    - consider adding an object
    - add distributions for lights?
        - not sure, since we are not necessarily doing the same thing
        - check the math for sampling
- cleanup sampling functions everywhere
    - probably removing sample_points/lines/triangles
    - cleanup sampling in ray tracing
    - make lights with single shapes in trace
- add radius in offsetting rays
- sample background to sum all environments
- envmap sampling
    - simplest case, pick pixels as 1D distribution
- Simple denoiser
    - joint bilateral denoiser
    - non-local means denoiser
    - denoiser from paris
- add builtin envmap to trace

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

## One shape

- change shape to use constant radius, fixed color
- hairball scene needs splitting for now
- update list marks shape buffers
- scene with name
- tesselation takes tags
- OpenGL with multiple index buffers
    - OpenGL updates: rebuild all buffers or detect if same size
    - selection carries shape ids
- BVH with multiple primitives
- All functions take all primitives
- shape with type
- facet_shape and friends are not virtual in API

## Test scenes: simnplify shape generation

- substance-like shader ball
- squircle
- bent floor
- add cube based tesselation for cylinder, disk
- use welding instead of complex grids
- 0 roughness

## OpenGL/Trace

- investigate bump map on GPU
    - https://www.opengl.org/discussion_boards/showthread.php/162857-Computing-the-tangent-space-in-the-fragment-shader
    - http://jbit.net/~sparky/sfgrad_bump/mm_sfgrad_bump.pdf

## Scene Import

- remove faceted shading
- PBR in OBJ
    - http://exocortex.com/blog/extending_wavefront_mtl_to_support_pbr
    - Pr/map_Pr (roughness) // new
    - Pm/map_Pm (metallic) // new
    - Ps/map_Ps (sheen) // new
    - Pc (clearcoat thickness) // new
    - Pcr (clearcoat roughness) // new
    - Ke/map_Ke (emissive) // new
    - aniso (anisotropy) // new
    - anisor (anisotropy rotation) // new

## yView

- double sided in scene
    - check glTF write
    - extend OBJ
- bug in light detection
    - check mcguire car
- cutout not working on yview
    - consider cutout by default
- add print scene info button to yview/ytrace
- add view cam settings

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


## Deployment

- shorter doc formatting
- postponed: Doxygen and Sphynx
- postponed: consider amalgamation for yocto and ext

## Math

- consider types without constructors
- consider removing const refs
- make make_basis
- span
- make_vec
- generic transform with make_vec and project_homogeneous
- frame inverse with flag
- make stronger the assumption on the use of frames
    - frame inversion
    - documentation
- check random shuffle
- check random number generation for float/double
- check rotation and decompoaition of rotations
   - see euclideanspace.com

## Scene

- add material to env
- remove node children: use stable sort
- share texture info accross GPU/tracer/scene
- make texture info more complete with mirroring and mipmapping
- update convert functions to new api (?)
- cleanup tesselation in shape
    - remove tesselate once
    - tesselation uses only internal levels

## Image

- maybe: make image a simple structure
    - get_pixel, make_image
- remove constructors and accessors from vec/mat/frame

## Low-level code

- serialization with visitor?
    - decide if exposing json is reasonable (no, move to C++17 any when possible)
      for now this is just a matter of compilation time
      later it is best to use a variant type
- serialization to JSON as the main reflection model

## Ui

- update imgui and imgui opengl3
    - consider gl3w since ImGui integrates it and it is small
- consider not using sliders
- cleanup scene widgets
- add angle semantic
- add rotation
- add frame editing with decomposition
- add labels 2,3,4
- consider disabling visitors and using JSON serialization with schema-like things

## yScnProc

- general fixup
- print info
