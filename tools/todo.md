# Notes on future improvements of Yocto/GL

This file contains notes on future improvements of Yocto.
Please consider this to be just development notes and not any real planning.

## New scene

- add material to env
- remove node children
    - use stable sort
    - add local frame
- nodes point to shapes
    - instances are created on the fly
    - update hierarchy creates instances
- procedural instances are deleted
- instances are optional
    - rename instances for now
    - make drawing code work with only shapes
    - instances are added during node updates
    - remove instances from being central
    - instances are added during update_hierarchy()
    - loaders have no option to add or remove hierarchy
- remove instances
    - trace should not need instances
        - bvh should handle instances in a different way,
            - maybe just returning frames instead of instance ids
            - or return sid, iid as before
    - move frame to shape_group
    - use nodes everywhere instances where needed
- move to tagged shape

## Trace

- environment map with material
- remove instances from tracer
    - handle environments as missing shape or as a special shape
    - add special shape types: inf sphere (env) distant points
    - handle light as frame + shape (none for env) + ke + ke_txt
- envlight parametrization
    - bad weight for envmap
    - bad envmap rendering
- trace options
    - no mis
    - no env lights
- fast distribution sampling
- fresnel in brdf
    - rescale fresnel with roughness
    - fresnel in coefficients
    - fresnel in weights
- add shape methods
    - surface/curve/point
    - quads/points/triangles etc
- try hard to eliminate deltas
    - I do not think they actually work right
    - put stringent epsilons
- path tracer with mis
    - possible bug in light weight

## Refactor

- shape with type
- make make_basis
- facet_shape and friends are not virtual in API
- span
- material in shape_group
- cleanup
    - make_vec
    - generic transform with make_vec and project_homogeneous

## Internal

- move away from special functions in BVH
    - always use sort
    - provide a sort buffer

## Scene Import

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
- filmic tonemapping from tungsten
- obj
    - decide whether to put it in eval_norm (probably not)
    - insert facet shape call after loading?
    - tesselate may create smoothed vertices or during vbo creation
    - trace respect smoothing
- trace
    - bug in reflection delta
    - remove wo from point
    - double sided rendering in the brdf and emission
    - envmap point is just a point far away and normal pointing in
    - SAH based build
- double sided option in render params
    - this forces double sided over the material settings
- doule sided in scene <-> glTF
- bug in light detection
    - check mcguire car
- cutout not working on yview
- consider double sided by default
    - check pbrt
- consider cutout by default
- add cutout to trace

- add print scene info button to yview/ytrace
- add view cam settings
- add bbox to trace
- add builtin envmap to trace
- better default material
- better eyelight
- add prefiltered look to trace/view ?

## BVH

- SAH based build
- simplify build code

## Deployment

- shorter doc formatting
- postponed: Doxygen and Sphynx
- postponed: consider amalgamation for yocto and ext

## Math

- consider constexpr
- consider types without constructors
- consider removing const refs

## Scenes

- setup scene repo
- begin importing scenes
- make 4 scene variants
    - original
    - fixed
    - converted OBJ
    - converted GLTF
- consider putting OBJ extensions into its own files?

## Postponed: Shade uses render buffers

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
