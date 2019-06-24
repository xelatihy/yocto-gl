# Notes on material modeling in Yocto/Scene and Yocto/Trace

There are two models we can adopt for materials in Yocto/GL: monolitic generic 
material or a use an enum to specify a material type. In this document we 
review both options by listing how to represent common materials.

## Common Materials

Here is a list of common materials we want to render:

- matte: double-sided diffuse
- plastic: double-sided fresnel blend between diffuse and specular
- metal: double-sided fresnel specular
- car-paint: double-sided fresnel metal with double-sided fresnel coat
- leaves: double-sided diffuse with inverted diffuse (2 lobes), possibly with fresnel specular coat
- thin-glass: double-sided fresnel blend between reflection and transmission
- glass: single-sided fresnel blend between reflection and transmission; should do refraction and tir
- volume-glass: same as glass, but also activate volume scattering; we should probably always do this
- subsurface-scattering: like glass but has volume albedo

## Material with Types

The simplest manner to represent the above materials is to use a material with 
separate types. Here is a list of parameters needed. We use the term 'base' to 
indicate the base color of the object and 'specular` to indicate the dieletric
coating.

- matte: base
- plastic: base, specular, roughness
- metal: base, roughness
- car-paint: base, specular, roughness
- leaves: base, subsurface, specular, roughness
- thin-glass: base (glass transmission color), specular, roughness
- glass: base (glass transmission color), specular, roughness
- volume-glass: base (color of transmissionm at 0.01m), specular, roughness
- subsurface-scattering: base (color of transmissionm at 0.01m), subsurface 
    (volume scattering), specular, roughness, phase g, volume depth

UI/File/Scene Parameters:

- base color
- specular color
- roughness
- subsurface color (phase g and volume depth, maybe)
- specular ior (alternative considered below)

Concers:

- fresnel is coubled with specular intensity
    - pros: cannot make incorrect setups
    - cons: surfaces may be too reflecting at grazing angles
    - alternative: add a specular ior only for specular, but not metals
  
Pros:

- simple to understand
- smaller number of parameters to setup
- implementation is decoupled; e.g. so if you fix a bug in glass, you are not   
    breaking plastics; in practice this was a major issue in getting volumes 
    to work (and it still is), and it also affect refracting glass
- with the exception of leaves and subsurface, all materials send out only two
    directions
- might be easier to extend,  but we do not know

Cons:

- more code path to write and test
- possibly harder to translate to/from OBJ and glTF formats

## Monolitic Material

As an alternative, chosen in film production, is to use a material with a bunch
of parameters that can be setup to produce the wanted results.

- matte: diffuse, thin
- plastic: diffuse, specular, roughness, thin
- metal: specular, roughness, thin; use edge_color to derive complex ior
- car-paint: specular, coat, roughness, thin
- leaves: diffuse, subsurface, specular, roughness, thin
- thin-glass: transmission (glass transmission color), specular, roughness
- glass: transmission (glass transmission color), specular, roughness
- volume-glass: transmission (surface), voltransmission (color of transmission
    at 0.01m), specular, roughness
- subsurface-scattering: transmission (surface), voltransmission (color of
    transmissionm at 0.01m), subsurface (volume scattering), specular,
    roughness, phase g, volume depth

Concers:

- everything might work is we flip the normal only for thin surfaces
  - cons: not implemented yet, so likely big stuff will break somewhere
  - cons: for sampling/pdf code it is hard to understand what this flip do
  - bottom line: we can try this and it may or may not work, but will certainly
    require a lot of coding to get there and retest everything
- fresnel is coubled with specular intensity
  - pros: cannot make incorrect setups
  - cons: surfaces may be too reflecting at grazing angles
- alternative: add a specular ior only for specular, but not metals
  - pros: above
  - cons: tir does not make sense in this case if specular != (1, 1, 1)
  
Pros:

- simpler i/o, ui and system code since no switches are necessary
- relatively strightforward to create new behaviours by adding additional lobes

Cons:

- lobes are coupled by fresnel and this makes it hard to make sure all
  cases are propertly handled in sampling and evaluation; pbrt special cases
  code for this going outside of its sum-of-lobes approach and switch to 
  a lobe "tree" so that fresnel mixtures are handled separately
- may be hard to properly set all parameters since they are couples, e.g.
  transmission, voltransmission, volscattering need to ve carefully setuo
- in most cases, only a few lobes are active, but the rendered needs to 
  handle them all since they could all be active
- to get all possible behaviours, either "smartly" configure material 
  (essentially reverse engineer the renderer) or add more parameters
  (production shaders have tons to configure all possible cases)


## Other issues

- pain points now
  - tir does not work and it is hard to make it work in the current framework
  - we need to keep both thin and non-thin materials, delta and non-delta
  - coupling between lobes makes it so freaking hard to change this code
    - i had multiple times fix one scene to see another break in a non-obvious
        manner; this is hallmark of unmaintainable code
- deltas
  - so far, no good way to handle deltas in shader without blurring, 
    especially for glass; maybe this could be a small paper to write MIS so
    that it always work with deltas
