# Yocto/Scene: Scene representation

Yocto/Scene define a simple scene representation, and related utilities,
used in the Yocto/GL path tracer and for scene IO.
Yocto/Scene is implemented in `yocto_scene.h` and `yocto_scene.cpp`.

## Scene representation

Scenes are stored in `scene_model` structs and are comprised of arrays of objects.
Scenes are comprised of camera, instances, shapes, materials, textures
and environments, and stored as arrays named as before.
Animation is not currently supported.
The scene representation is geared toward modeling physically-based environments.
In Yocto/Scene, lights are not explicitly defined, but implicitly comprised of
instances with emissive materials and environment maps.
All scenes and objects properties are accessible directly.

All scenes objects may have names that are used in IO. If names are defined,
that have to be unique. If not, names are automatically generated. Names are
stored separately from objects, for performance reasons. So for each object
array, Yocto/Scene stores a corresponding names array. For examples,
cameras as stored as `cameras` and their names are stored as `camera_names`.

Cameras, instances and environments have coordinate frames to define the local
to world transformation. Frames are presented as affine 3x4 matrices and are
intended to be rigid transforms, although most scene processing support frames
with scaling.

**Cameras**, represented by `scene_camera`, are based on a simple lens model.
Cameras coordinate systems are defined by their frame.
Cameras projections are described in photographic terms. In particular,
we specify film size (35mm by default), film aspect ration,
the lens' focal length, the focus distance and the lens aperture.
All values are in meters. We support both perspective and orthographic cameras,
but prefer the former.

Common aspect ratios used in video and still photography are
3:2 on 35 mm (0.036 x 0.024),
16:9 on 35 mm (0.036 x 0.02025 or 0.04267 x 0.024),
2.35:1 on 35 mm (0.036 x 0.01532 or 0.05640 x 0.024),
2.39:1 on 35 mm (0.036 x 0.01506 or 0.05736 x 0.024),
2.40:1 on 35 mm (0.036 x 0.015 or 0.05760 x 0.024).
To compute good apertures, one can use the F-stop number from photography
and set the aperture to focal length over f-stop.

**Textures**, represented as `scene_texture`, contains either 8-bit LDR or
32-bit float HDR images with four channels. Textures can be encoded in either
a linear color space or as sRGBs, depending on an internal flag. The use of
float vs byte is just a memory saving feature.

**Materials**, represented as `scene_material`, are defined by a material type
and a few parameters, common to all materials. In particular, we support the
following materials:

- `matte`, for materials like concrete or stucco, implemented as a lambertian bsdf;
- `glossy`, for materials like plastic or painted wood, implemented as the sum
  of a lambertian and a microfacet dielectric lobe;
- `metallic`, for materials like metals, implemented as either a delta or
  microfacet brdf lobe;
- `transparent`, for materials for thin glass, implemented as a delta or
  microfacet transmission bsdf;
- `refractive`, for materials for glass or water, implemented as a delta or
  microfacet refraction bsdf; also support homogenous volume scattering;
- `subsurface`, for materials for skin, implemented as a microfacet refraction
  bsdf with homogenous volume scattering - for no this is like `refractive`;
- `volume`, for materials like homogeneous smoke or fog, implemented as the lack
  of a surface interface but with volumetric scattering.
- `gltfpbr`, for materials that range from glossy to metallic, implemented as
  the sum of a lambertian and a microfacet dielectric lobe;
  this is a compatibility material for loading and saving Khronos glTF data.

All materials can specify a diffuse surface emission `emission` with HDR values
that represent emitted radiance.

Surface scattering is controlled by specifying the main surface color `color`,
that represent the surface albedo, the surface roughness `roughness` and
the index of refraction `ior`. The physical meaning of each parameter depends
on the material type. By default surfaces are fully opaque, but
can defined a `opacity` parameter and texture to define the surface coverage.

Materials like `refractive`, `subsurface` and `volume` may also specify
volumetric properties. In these cases, the `color` parameter controls the volume density,
while the `scattering` also define volumetric scattering properties by setting a
`transmission` parameter controls the homogenous volume scattering.

All parameters can modulated by a corresponding textures, if present.

**Shapes**, represented by `scene_shape`, are indexed meshes of elements.
Shapes can contain only one type of element, either
points, lines, triangles or quads. Shape elements are parametrized as in
[Yocto/Geometry](yocto_geometry.md).
Vertex properties are defined as separate arrays and include
positions, normals, texture coords, colors, radius and tangent spaces.
Additionally, Yocto/Scene supports face-varying primitives, as `scene_fvshape`,
where each vertex data has its own topology.

**Instances**, represented as `scene_instance`, place shapes in the scene by
defining their coordinate frame, a shape index and a material index.
Through the use of instancing, Yocto/Scen scales well to large environments
without introducing more complex mechanisms.

**Environments**, represented as `scene_environment`, store the background
illumination as a scene. Environments have a frame, to rotate illumination,
an emission term and an optional emission texture.
The emission texture is an HDR environment map stored in a LatLon
parametrization.

**Subdivs**, represented as `scene_subdiv`, support tesselation and displacement
mapping. Subdivs are represented as facee-varying shapes.
Subdivs specify a level of subdivision and can be subdivide elements
either linearly or using Catmull-Clark subdivision. Subdivs also support
displacement by specifying both a displacement texture and a displacement amount.
Differently from most systems, in Yocto/Scene displacement is specified
in the shape and not the material. Subdivs only support tesselation to shapes,
but do not directly support additional evaluation of properties.
Subdivs specified to the shape index to which they are sub divided into,
and provide tesselation support as discussed later.

## Scene Creation

Objects are added to the scene by directly adding elements to the corresponding
arrays. References to elements are expressed as indices to the
corresponding arrays. For each element type, properties can be set directly.

For cameras, you should set the camera frame, the camera view,
via lens, aspect and film, and optionally the camera aperture and focus.

```cpp
auto scene = scene_model{};                  // create a scene
auto& camera = scene.cameras.emplace_back(); // create a camera
camera.frame = identity3x4f;                 // set frame to identity
camera.lens = 0.050;                         // set as 50mm lens
camera.aspect = 1.5;                         // set 3:2 aspect ratio
camera.film = 0.036;                         // set the film as 35mm
camera.aperture = 0.01;                      // set 10mm aperture
camera.focus = 10;                           // set the focus at 10m
```

For instances, you should set the instance frame, shape and material.

```cpp
auto scene = scene_model{};                         // create a scene
auto& instance = scene.instances.emplace_back();    // create an instance
instance.frame = identity3x4f;                      // set frame to identity
auto& shape = scene.shapes.emplace_back();
instance.shape = (int)scene.shapes.size()-1;        // set shape index
auto material = add_material(scene);
instance.material = (int)scene.materials.size()-1;  // set material index
```

For textures, set the size, the color space, and _either_ the hdr or ldr pixels.

```cpp
auto scene = scene_model{};                   // create a scene
auto texture = scene.textures.emplace_back(); // create a texture
texture.width = ...;  texture.height = ...;   // set size
texture.linear = ...;                         // set color space
if (...) {
  texture.pixelsf = vector<vec4f>{...};       // set float pixels
} else {
  texture.pixelsb = vector<vec4b>{...};       // set byte pixels
}
```

For materials, we need to specify the material type and color at the minimum.
We can further control the appearance by changing surface roughness, index of
refraction and volumetric properties, when appropriate. Here are some examples.

```cpp
auto scene = scene_model{};                   // create a scene
auto& matte = scene.materials.emplace_back(); // create a material
matte.type = scene_material_type::matte;      // with matte appearance
matte.color = {1,0.5,0.5};                    // with base color and
matte.color_tex = texture_id;                 // textured albedo
auto& glossy =  scene.materials.emplace_back(); // create a material
glossy.type = scene_material_type::glossy;    // with glossy appearance
glossy.color = {0.5,1,0.5};                   // with constant color
glossyv.roughness = 0.1;                      // base roughness and a
glossy.roughness_tex = add_texture(scene);    // roughness texture
auto& metallic =  scene.materials.emplace_back(); // create a material
glossy.type = scene_material_type::glossy;    // with metallic appearance
metal.color = {0.5,0.5,1};                    // constant color
metal.roughness = 0.1;                        // constant roughness
auto& tglass = scene.materials.emplace_back(); // create a material
tglass.type = scene_material_type::transparent;// with a transparent appearance
tglass.color = {1,1,1};                       // with constant color
auto& glass = scene.materials.emplace_back(); // create a material
glass.type = scene_material_type::transparent;// with a refractive appearance
glass.color = {1,0.9,0.9};                    // constant color
auto& subsurf = scene.materials.emplace_back();// create a material
subsurf.type = scene_material_type::subsurface;// with a refractive appearance
subsurf.color = {1,1,1};                      // that transmits all light
subsurf.scattering = {0.5,1,0.5};             // and as volumetric scattering
```

Lights are not explicit in Yocto/Scene but are specified by assigning emissive
materials.

```cpp
auto scene = scene_model{};                   // create a scene
auto& light = scene.materials.emplace_back(); // create a material
light.color = {0,0,0};                        // that does not reflect light
light.emission = {10,10,10};                  // but emits it instead
```

For shapes, you should set the shape elements, i.e. point, limes, triangles
or quads, and the vertex properties, i.e. positions, normals, texture
coordinates, colors and radia. Shapes support only one element type.

```cpp
auto scene = scene_model{};                   // create a scene
auto& shape = scene.shapes.emplace_back();    // create a shape
shape.triangles = vector<vec3i>{...};         // set triangle indices
shape.positions = vector<vec3f>{...};         // set positions
shape.normals = vector<vec3f>{...};           // set normals
shape.texcoords = vector<vec2f>{...};         // set texture coordinates
```

We also support subdivision surfaces, stored as face-varying shapes.
In this case, set the quads for positions, normals and texture coordinates.
Also set the subdivision level, and whether to use Catmull-Clark or linear
subdivision. Finally, displacement can also be applied by setting a displacement
scale and texture.

```cpp
auto scene = scene_model{};                     // create a scene
auto& subdiv = scene.subdivs.emplace_back();    // create a subdiv
subdiv.quadspos = vector<vec4i>{...};           // set face-varying indices
subdiv.quadstexcoord = vector<vec4i>{...};      // for positions and textures
subdiv.positions = vector<vec3f>{...};          // set positions
subdiv.texcoords = vector<vec2f>{...};          // set texture coordinates
subdiv.subdivisions = 2;                        // set subdivision level
subdiv.catmullclark = true;                     // set Catmull-Clark subdivision
subdiv.displacement = 1;                        // set displacement scale
subdiv.displacement_tex = texture_id;           // and displacement map
```

For environments, set the frame, emission and optionally the emission texture.

```cpp
auto scene = scene_model{};                     // create a scene
auto& environment = scene.environments.emplace_back(); // create an environment
environment.frame = identity3x4f;               // set identity transform
environment.emission = {1,1,1};                 // set emission scale
environment.emission_tex = texture_id;          // add emission texture
```

We also support face-varying shapes, albeit not stored in the scene
(see subdivs above). In this case, set the quads for positions,
normals and texture coordinates.

```cpp
auto& shape = scene_fvshape{};              // create a shape
shape.quadspos = vector<vec4i>{...};        // set face-varying indices
shape.quadstexcoord = vector<vec4i>{...};   // for positions and textures
shape.positions = vector<vec3f>{...};       // set positions
shape.texcoords = vector<vec2f>{...};       // set texture coordinates
```

## Scene tesselation

Scene specify geometry as instance shapes. For convenience we also support
subdivision surfaces, that need to be tessellated into indexed meshes before
use with `tesselate_sudivs(scene)`.

```cpp
auto scene =  scene_model{...};          // create a complete scene
void tesselate_subdivs(scene);           // tesselate subdivs
```

## Evaluation of scene properties

Yocto/Scene defines several function to evaluate scene properties.
Use `compute_bounds(scene)` to compute the scene bounding boxes.
Use `get_camera(scene, name)` to get a camera by name or the default camera
is the name is not given. Use `eval_camera(camera, image_uv, lens_uv)`
to get a camera ray from the normalized image coordinates `image_uv` and
lens coordinates `lens_uv`.

```cpp
auto scene = new trace_scene{...};             // create a complete scene
auto camera = get_camera(scene);               // get default camera
auto ray = eval_camera(camera,{0.5,0.5},{0,0});// get ray though image center
```

Use `texture_size(texture)` to get the texture resolution, and
`eval_texture(texture, uv)` to evaluate the texture at specific uvs.
Textures evaluation returns a color in linear color space, regardless of
the texture representation.

```cpp
auto scene = new trace_scene{...};           // create a complete scene
auto texture = scene.textures.front();        // get first texture
auto col = eval_texture(texture,{0.5,0.5});    // eval texture
```

Use `eval_material(material, texcoord)` to evaluate material textures and
combine them with parameter values. The function returns a
`scene_material_sample` that has the same parameters of a material but no
textures defined.

```cpp
auto scene = new trace_scene{...};             // create a complete scene
auto material = scene.materials.front();        // get first material
auto mat = eval_material(material,{0.5,0.5});    // eval material
```

Several functions are defined to evaluate the geometric and material
properties of points on instances, indicated by the shape element id
and, when needed, the shape element barycentric coordinates.
Use `eval_position(...)` to evaluate the point position,
`eval_normal(...)` to evaluate the interpolate point normal,
`eval_texcoord(...)` to evaluate the point texture coordinates,
`eval_element_normal(...)` to evaluate the point geometric normal, and
`eval_color(...)` to evaluate the interpolate point color.
Use `eval_material(...)` as a convenience function to evaluate material
properties of instance points.

```cpp
auto scene = new trace_scene{...};             // create a complete scene
auto instance = scene.instances.front();      // get first instance
auto eid = 0; auto euv = vec3f{0.5,0.5};       // element id and uvs
auto pos  = eval_position(instance, eid, euv); // eval point position
auto norm = eval_normal(instance, eid, euv);   // eval point normal
auto st   = eval_texcoord(instance, eid, euv); // eval point texture coords
auto col  = eval_color(instance, eid, euv);    // eval point color
auto gn   = eval_element_normal(instance, eid, euv); // eval geometric normal
auto mat  = eval_material(instance, eid, euv); // eval point material
```

Use `eval_environment(environment, direction)` to evaluate an environment
map emission along a specific direction `direction`. Use
`eval_environment(scene, direction)` to accumulate the lighting for all
environment maps.

```cpp
auto scene = new trace_scene{...};               // create a complete scene
auto enva = eval_environment(scene, dir);        // eval all environments
auto environment = scene.environments.front();  // get first environment
auto envi = eval_environment(environment, dir);  // eval environment
```

## Serialization formats

Yocto/SceneIO supports loading and saving to Ply, Obj, Pbrt, glTF,
and a custom Json format. For the standard formats, loading is best effort,
since scene data is transformed from the formats' scene models to the
Yocto/SceneIO model.

The custom Json format is a serialization of the internal properties for
most scene objects, with a few conventions taken for extensibility.
Scene's object arrays are represented as dictionaries in Json with the
objects' names used as keys. This ensure proper reference semantic and
allows for more extensibility in the future, but it also means that
object order is not preserved during serialization.

Cameras, materials, instances and environments are represented directly
in the Json scene, while textures and shapes are serialized using
standard image and geometry formats. By convention, scenes are
stored as a single Json format for the scene structure. Textures are
stored in the `textures` directory with the name of the texture as filename,
while the extension is determined by checking th available files.
Shapes are stored in the `shapes` directory with name of the shape as filename,
while the extension is determined by checking th available files.

## Loading and saving scenes

Scenes are loaded with `load_scene(filename, scene, error, progress)` and
saved with `save_scene(filename, scene, error, progress)`.
Both loading and saving take a filename, a scene pointer and return
whether or not the scene was loaded successfully.
In the case of an error, the IO functions set the `error` string with a
message suitable for displaying to a user.
The functions take a progress callback as an optional parameter,
that is called as scene loading progresses.

```cpp
auto scene = new sceneio_scene{};                    // scene
auto progress = [](const string& message,          // progress callback
                   int current, int total) {
  print_info(message, current, total);
};
auto error = string{};                             // error buffer
if(!load_scene(filename, scene, error, progress))  // load scene
  print_error(error);
if(!save_scene(filename, scene, error, progress))  // save scene
  print_error(error);
```

## Scene stats and validation

Yocto/SceneIO has functions to compute scene stats and provide validation of
scene data. Use `scene_stats(scene)` to get scene stats and
`scene_validation(scene)` to validate scene objects.

```cpp
auto scene = new sceneio_scene{...};          // create a complete scene
auto stats = scene_stats(scene);            // get stats
for(auto stat : stats) print_info(stat);    // print stats
auto errors = validate_stats(scene);        // get validation errors
for(auto error : errors) print_error(error);// print error
```

## Example scenes

Yocto/SceneIO has a function to create a simple Cornell Box scene for testing.
There are plans to increase support for more test scenes in the future.

```cpp
auto scene = new sceneio_scene{...};          // create a complete scene
make_cornellbox(scene);                     // make cornell box
```

## Evaluation of scene properties

Yocto/SceneIO defines several function to evaluate scene properties.
Use `compute_bounds(scene)` to compute the scene bounding boxes.
Use `get_camera(scene, name)` to get a camera by name or the default camera
is the name is not given. Use `eval_camera(camera, image_uv, lens_uv)`
to get a camera ray from the normalized image coordinates `image_uv` and
lens coordinates `lens_uv`.

```cpp
auto scene = new sceneio_scene{...};             // create a complete scene
auto camera = get_camera(scene);               // get default camera
auto ray = eval_camera(camera,{0.5,0.5},{0,0});// get ray though image center
```

Use `texture_size(texture)` to get the texture resolution, and
`eval_texture(texture, uv)` to evaluate the texture at specific uvs.
Textures evaluation returns a color in linear color space, regardless of
the texture representation.

```cpp
auto scene = new sceneio_scene{...};           // create a complete scene
auto texture = scene.textures.front();        // get first texture
auto col = eval_texture(texture,{0.5,0.5});    // eval texture
```

Use `eval_material(material, texcoord)` to evaluate material textures and
combine them with parameter values. The function returns a
`scene_material_sample` that has the same parameters of a material but no
textures defined.

```cpp
auto scene = new sceneio_scene{...};             // create a complete scene
auto material = scene.materials.front();        // get first material
auto mat = eval_material(material,{0.5,0.5});    // eval material
```

Several functions are defined to evaluate the geometric and material
properties of points on instances, indicated by the shape element id
and, when needed, the shape element barycentric coordinates.
Use `eval_position(...)` to evaluate the point position,
`eval_normal(...)` to evaluate the interpolate point normal,
`eval_texcoord(...)` to evaluate the point texture coordinates,
`eval_element_normal(...)` to evaluate the point geometric normal, and
`eval_color(...)` to evaluate the interpolate point color.
Use `eval_material(...)` as a convenience function to evaluate material
properties of instance points.

```cpp
auto scene = new sceneio_scene{...};             // create a complete scene
auto instance = scene.instances.front();      // get first instance
auto eid = 0; auto euv = vec3f{0.5,0.5};       // element id and uvs
auto pos  = eval_position(instance, eid, euv); // eval point position
auto norm = eval_normal(instance, eid, euv);   // eval point normal
auto st   = eval_texcoord(instance, eid, euv); // eval point texture coords
auto col  = eval_color(instance, eid, euv);    // eval point color
auto gn   = eval_element_normal(instance, eid, euv); // eval geometric normal
auto mat  = eval_material(instance, eid, euv); // eval point material
```

Use `eval_environment(environment, direction)` to evaluate an environment
map emission along a specific direction `direction`. Use
`eval_environment(scene, direction)` to accumulate the lighting for all
environment maps.

```cpp
auto scene = new sceneio_scene{...};               // create a complete scene
auto enva = eval_environment(scene, dir);        // eval all environments
auto environment = scene.environments.front();  // get first environment
auto envi = eval_environment(environment, dir);  // eval environment
```
