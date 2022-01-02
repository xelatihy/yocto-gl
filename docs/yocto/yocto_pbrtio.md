# Yocto/PbrtIO: Serialization for Pbrt models

Yocto/ModelIO is a collection of utilities for loading and saving scenes
in the Pbrt format.
Yocto/PbrtIO is implemented in `yocto_pbrtio.h` and `yocto_pbrtio.cpp`,
and depends on `fast_float.h`.

## Pbrt models

The Pbrt file format is a scene representation suitable for realistic rendering
and implemented by the Pbrt renderer.
To use this library is helpful to understand the basic of the Pbrt
file format for example from the
[Pbrt format documentatioon](https://www.pbrt.org/fileformat-v3.html).

The Pbrt file format is an extensible file format for a plugin-based system.
Representing the format directly allows for best fidelity but pushed the burden
of interpreting standard plugins to the use. Yocto/ModelIO takes a different
approach and translates camera, shapes, materials, textures and lights from
Pbrt plugins to a common representation that presents users a simpler and
more uniform scene representation.

Yocto/ModelIO represents Pbrt data with the `pbrt_scene` struct.
Pbrt models are defined as collections of cameras, instanced shapes, materials,
texture and environments. Pbrt cameras are translate into a thin-len
approximations. Pbrt materials are translated to a material representation
similar to the Disney BSDF. Pbrt textures are either interpreted ion the fly or
defined by a image file. Pbrt area lights are translated to either emissive
materials and environments.

The Pbrt model is defined as an array of objects of the types defined above.
Pbrt objects are pointers owned by the main `pbrt_scene`.
Objects properties can be read and written directly from the model data,
and are documented in the header file for now.
Yocto/ModelIO does not currently provide functions to read and write Pbrt
shapes with a simpler interface than accessing data directly.

In general, Pbrt support is still experimental even if the library can
parse most Pbrt files. The objects properties documentations are for now
stored in the header file.

```cpp
auto pbrt = new pbrt_scene{...};            // obj model buffer
for(auto shape : pbrt.shapes)              // access shapes
  print_info(shape.name);                  // access shape properties
for(auto material : pbrt.material)         // access materials
  print_info(material.diffuse);            // access material properties
for(auto material : pbrt.material)         // access materials
  print_info(material.color_tex);          // access material textures
for(auto camera : pbrt.cameras)            // access cameras [extension]
  print_info(camera.frame);                // access camera properties
for(auto environment : pbrt.environments)  // access environments [extension]
  print_info(environment.emission);        // access environment properties
```

Use `ok = load_pbrt(filename, pbrt, error)` to load Pbrt 
files and `ok = save_pbrt(filename, pbrt, error)` to save them.  
Both loading and saving take a filename, a reference to a model, 
and an error string, and returns a boolean that indicated whether the 
operation is successful.

```cpp
auto pbrt = pbrt_model{};              // pbrt model
auto error = string{};                 // error string
if (!load_pbrt(filename, pbrt, error)) // load pbrt
  handle_error(error); 
if (!save_pbrt(filename, pbrt, error)) // save pbrt
  handle_error(error); 
```
