# Notes on future improvements of Yocto/GL

This file contains notes on future improvements of Yocto.
Please consider this to be just development notes and not any real planning.

## Main features

- better tests
- better material rendering
- better rendering
- interactive procedural shapes
- prepare for research on procedural components

## Port scenes

- yocto
    - load pfm
    - pixel clamp
    - pixel filter
    - fresnel scale with roughness (should be option)
    - fresnel scale transmission
    - tone enhancement
    - fresnel sampling
    - better bvh?
    - highlights are too soft in bitterli scenes
- bitterli
    - render exclude list
    - bathroom
    - bathroom2
    - bedroom
    - classroom
        - light
    - coffee
    - dining-room
        - light
    - dragon
        - light
    - house
        - light
    - kitchen
    - lamp
        - light
    - living-room
        - check render
    - living-room-2
        - check render
    - living-room-3
        - check render
    - spaceship
        - check render
    - staircase
    - staircase2
    - teapot
- bitterli
    - hair
    - veach-bidir
        - bad materials
        - too high energy
    - kitchen
        - bump map
    - dining room
        - infinite sphere cap
    - lamp
        - bulb size
    - car
        - shading bug
        - enable kr
        - metallic paint
        - envmap
    - car2
        - shading bug
        - enable kr
        - metallic paint
    - check texture color
    - check tone mapping
    - metallic paint
        - car, car2
    - huge number of bounces
    - make list of scenes that should not be ported
    - flip double sided normals on at least large objects
- mcguire
    - list of scenes
    - some scenes have wrong transparency
    - cameras
    - lights
    - write list of issues
    - add bounding box print
- move ytrace batch size
    - use callback workflow
- pbrt
    - pber parser
- gltf exports

## Trace

- move ytrace batch size
- roussian roulette on weight
- samplers
    - sobol
    - adaptive sampling ala tungsgen
- shape
    - faceted shading
- brdf
    - deltas without delta flag
    - check pbrt hair
- light sampling
    - possible bug in light weight
    - envmap sampling
    - path trace with explicit light sampling
    - check tungsten light smapling
    - eval_direct function
    - mis in params and not renderer?
- fresnel in brdf
    - fix kr
    - rescale fresnel with roughness
    - fresnel in coefficients
    - fresnel in weights

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
- shapes
    - squircle
    - bent floor
- add cube based tesselation for cylinder, disk
- use welding instead of complex grids
- 0 roughness
- transparent
- fix obj export
    - check shape names
    - save_obj()
        - skip group names if only one group
        - skip smoothing if all on

## Ui: clenanup scene widgets

- remove draw_value_widgets

## New scene

- add material to env
- remove node children
    - use stable sort
    - add local frame

## OpenGL/Trace

- optional post event on OSX, disable on Linux
- OpenGL new version 4.1
- use derivatives
    - for triangles, compute flat shading and triangle edges well
        - no need for explicit edges
        - http://www.aclockworkberry.com/shader-derivative-functions/
        - https://github.com/rreusser/glsl-solid-wireframe
- investigate bump map on GPU
    - https://www.opengl.org/discussion_boards/showthread.php/162857-Computing-the-tangent-space-in-the-fragment-shader
    - http://jbit.net/~sparky/sfgrad_bump/mm_sfgrad_bump.pdf

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
- move away from special functions in BVH?
    - always use sort
    - provide a sort buffer

## Deployment

- shorter doc formatting
- postponed: Doxygen and Sphynx
- postponed: consider amalgamation for yocto and ext

## Math

- consider constexpr
- consider types without constructors
- consider removing const refs
- make make_basis
- span
- make_vec
- generic transform with make_vec and project_homogeneous

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
