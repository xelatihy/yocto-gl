# Notes on future improvements of Yocto/GL

This file contains notes on future improvements of Yocto.
Please consider this to be just development notes and not any real planning.

## Main features

- better tests
- better material rendering
- better rendering

## planned cleanup

- more image operations to floats
    - make explicit the use of 1/3/4 channel images
    - bumps and dimples should be sphere
- single material representation converted on load
- convert textures on load
    - srgb -> float on load with correct convention

## glTF tests

- yview
    - set focus to bbox center
- yitrace
    - envmap rendering
- update normals/tangents
- yview
    - 2CylinderEngine: near plane
    - AlphaBlendModeTest: not supported
    - AnimatedMorphCube/AnimatedMorphSphere: not supported
    - AnimatedTriangle: crash
    - Avocado: bad specular?
    - BarracubeFish: flip normal map? specular?
        - occlusion not supported
    - BoomBox: specular?
    - BoxVertexColor: check alpha?
    - Buggy: transforms
    - DamagedHelmet: no colors
    - Duck: not visible 
    - MetalRoughSphere: check with envmap
    - NormalTangentMirrorTest: bad normal and bad sideness
    - NormalTangentTestL check colors
    - SimpleMeshes: crash
    - SimpleMorph: crash
    - SimpleSparseAccessor: crash
    - Suzanne: bad materials
    - Triangle: crash
    - TextureSettingTest: bad texture and double-sided
    - VC: bad transforms

## Trace

- pass gltf model
- bump/normal/displacement mapping
    - convert bump maps to normal maps internally
    - eval tangent space per triangle
    - compute derivaties for textures
- hair example
- subdiv example
- put back double sided
- refraction
    - rough surface refraction
- redirect log to imgui
- points/lines projections
- add radius in offsetting rays
- bvh with opacity
- simple denoiser
    - joint bilateral denoiser
    - non-local means denoiser
    - denoiser from paris

## Simpler shapes

- sampling shapes
    - remove sample_points, sample_lines, etc...
    - support uniform sampling in sample_shape

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
