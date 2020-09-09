# Yocto/SceneIO: Scene serialization

Yocto/SceneIO defines a simple scene representation, and related utilities,
mostly geared towards scene creation and serialization.
Yocto/SceneIO supports loading and saving scenes from Ply, Obj, Pbrt, glTF
and a custom Json format.
Yocto/SceneIO is implemented in `yocto_sceneio.h` and `yocto_sceneio.cpp`,
and depends on `cgltf.h`.

## Scene representation

Scenes are stored in `sceneio_scene` structs and are comprised of array
of objects whose memory is owned by the scene.
Scenes are comprised of camera, instances, shapes, materials, textures
and environments. Animation is not currently supported.
The scene representation is geared toward modeling physically-based environments
and and is similar to the [Yocto/Trace](yocto_trace.md) scene.
In Yocto/SceneIO, lights are not explicitly defined, but
implicitly comprised of instances with emissive materials and environment maps.
All scenes and objects properties are accessible directly.
Some scene data, like lights and ray-acceleration structures,
are computed when necessary and should not be accessed directly.

All scenes objects have names that have to be unique and
are used for both UI and serialization. Cameras, instances and environments
have coordinate frames to define the local to world transformation.
Frames are presented as affine 3x4 matrices and are intended to be
rigid transforms, although most scene processing support frames with
scaling.

**Cameras**, represented by `sceneio_camera`, are based on a simple lens model.
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

**Textures**, represented as `sceneio_texture` contain either 8-bit LDR or
32-bit float HDR images with four channels.
HDR images are encoded in linear color space, while LDR images
are encoded in sRGB.

**Materials** are modeled similarly to the
[Disney Principled BSDF](https://blog.selfshadow.com/publications/s2015-shading-course/#course_content) and the
[Autodesk Standard Surface](https://autodesk.github.io/standard-surface/).
Materials are defined using many parameters that control material emission,
surface scattering and homogeneous volumetric scattering.
Each material parameter has an associated texture, where texture values are
multiplied by parameter values.

Materials specify a diffuse surface emission `emission` with HDR values
that represent emitted radiance.

Surface scattering is modeled by defining the main surface
color `color`, that represent the surface albedo.
Specular reflection is modeled by default as a dielectric with index of
refraction `ior`, roughness `roughness` and is scaled by a specular
weight `specular`. By default, specular reflection is dielectric.
Materials can reflect light as metals by setting the `metallic` parameter.
In this case, a metallic specular reflection is defined by the surface color
interpreted as metallic reflectivity, and surface roughness defined by
the previous parameter. Surfaces are optionally covered by a coating layer
defined by a `coat` parameter. By default surfaces are fully opaque, but
can defined a `opacity` parameter and texture to define the surface coverage.

Materials also define volumetric scattering properties by setting a
`transmission` parameter that sets the transmitted surface color.
Surfaces are be default considered as thin sheet, but can be modeled as
isotropic volumes if the `thin` parameter is set to false. In that case,
the surface transmission controls the volumetric parameters by defining the
volume density, while the volume scattering albedo is defined by the
`scattering` property.

**Shapes** are represented as indexed meshes of elements using the
`sceneio_shape` type. Shapes can contain only one type of element, either
points, lines, triangles or quads. Shape elements are parametrized as in
[Yocto/Geometry](yocto_geometry.md).
Vertex properties are defined as separate arrays and include
positions, normals, texture coords, colors, radius and tangent spaces.
Additionally, Yocto/Scenes supports face-varying primitives where
each vertex data has its own topology.

Shapes supports tesselation and displacement mapping. Shape specify a
level of subdivision and can be subdivide elements either linearly
or using Catmull-Clark subdivision. Shapes also support displacement
by specifying both a displacement texture and a displacement amount.
Differently from most systems, in Yocto/SceneIO displacement is specified
in the shape and not the material, that only controls shading.
Subdivision and displacement are only specified in shapes, but not
taken into account when evaluating shape properties. For this to happen,
shapes have to be tessellated, as shown later.

Shapes are placed in the scene by defining shape **instances** that
take a coordinate frame, a shape pointer and a material pointer.
Instances are represented by the `sceneio_instance` type.
Instances are represented as `sceneio_instance` objects. Thought the
use of instancing Yocto/SceneIO scales well to large environments without
introducing more complex mechanisms.

Scenes might be lit by background illumination defined by **environments**,
represented by the `sceneio_environment` type. Environments have a frame,
to rotate illumination, an emission term and an optional emission texture.
The emission texture is an HDR environment map stored in a LatLon
parametrization.

## Scene Creation

Objects are added to the scene via `add_<object>(scene,name)` functions,
where `<object>` is the object type name. In these functions,
the name is optional and, if left blank, a unique name will be generated
automatically. For each object type, properties can be set directly.

For camera, you should set the camera frame, the camera view,
via lens, aspect and film, and optionally the camera aperture and focus.

```cpp
auto scene = new trace_scene{};       // create a scene
auto camera = add_camera(scene);      // create a camera named cam
camera->frame = identity3x4f;         // set frame to identity
camera->lens = 0.050;                 // set as 50mm lens
camera->aspect = 1.5;                 // set 3:2 aspect ratio
camera->film = 0.036;                 // set the film as 35mm
camera->aperture = 0.01;              // set 10mm aperture
camera->focus = 10;                   // set the focus at 10m
```

For instances, you should set the instance frame, shape and material.

```cpp
auto scene = new trace_scene{};       // create a scene
auto instance = add_instance(scene);  // create an instance named ist
instance->frame = identity3x4f;       // set frame to identity
auto shape = add_shape(scene);
instance->shape = shape;              // set shape pointer
auto material = add_material(scene);
instance->material = material;        // set material pointer
```

For textures, set _either_ the hdr or ldr image.

```cpp
auto scene = new trace_scene{};       // create a scene
auto texture = add_texture(scene);    // create a texture named tex
texture->hdr = image<vec4f>{...};     // set as a HDR texture
texture->ldr = ,image<vec4b>{...};    // set as a LDR texture
```

For materials, we adopt a Disney0like model that has many parameters,
but can render a large varierty of looks. Here are some examples.

```cpp
auto scene = new trace_scene{};               // create a scene
auto matte = add_texture(scene);              // create a matte material
matte->color = {1,0.5,0.5};                   // with baese color and
matte->color_tex = add_texture(scene);        // textured albedo
auto plastic = add_texture(scene);            // create a plastic material
plastic->color = {0.5,1,0.5};                 // with constant color
plastic->specular = 1;                        // constant specular
plastic->roughness = 0.1;                     // base roughness and a
plastic->roughness_tex = add_texture(scene);  // roughness texture
auto metal = add_texture(scene);              // create a metal material
metal->color = {0.5,0.5,1};                   // constant color
metal->metallic = 1;                          // constant metallic
metal->roughness = 0.1;                       // constant roughness
auto tglass = add_texture(scene);             // create a thin glass material
tglass->color = {1,1,1};                      // with constant color
tglass->specular = 1;                         // constant specular
tglass->transmission = 1;                     // constant transmission
auto glass = add_texture(scene);              // create a glass material
glass->color = {1,1,1};                       // constant color
glass->specular, = 1;                         // constant specular
glass->transmission = 1;                      // constant transmission
glass->thin = false;                          // non-volumetric material
auto subsurf = add_texture(scene);            // create a subsurface material
subsurf->color = {1,1,1};                     // constant color
subsurf->specular = 1;                        // constant specular
subsurf->transmission = 1;                    // constant transmission
subsurf->thin = false;                        // volumetric material
subsurf->scattering = {0.5,1,0.5};            // volumetric scattering
```

For shapes, you should set the shape elements, i.e. point, limes, triangles
or quads, and the vertex properties, i.e. positions, normals, texture
coordiantes, colors and radia. Shapes support only one element type.

```cpp
auto scene = new trace_scene{};             // create a scene
auto shape = add_shape(scene);              // create a shape named shp
shape->triangles = vector<vec3i>{...};      // set triangle indices
shape->positions = vector<vec3f>{...};      // set positions
shape->normals = vector<vec3f>{...};        // set normals
shape->texcoords = vector<vec2f>{...};      // set texture coordinates
```

Shapes can also be face-varying. In this case, set the quads for positions,
normals and texture coordinates. This is helpful when using subdivision
surfaces, which are specified by settings the subdivision level, and whether
to use Catmull-Clark or linear subdivision. Finally, displacement can also
be applied by setting a displacement scale and texture.

```cpp
auto scene = new trace_scene{};             // create a scene
auto shape = add_shape(scene);              // create a shape named shp
shape->quadspos = vector<vec4i>{...};       // set face-varying indices
shape->quadstexcoord = vector<vec4i>{...};  // for positions and textures
shape->positions = vector<vec3f>{...};      // set positions
shape->texcoords = vector<vec2f>{...};      // set texture coordinates
shape>subdivisions = 2;                     // set subdivision level
shape->catmullclark = true;                 // set Catmull-Clark subdivision
shape->displacement = 1;                    // set displacement scale
shape->displacement_tex = tex;              // and displacement map
```

For environments, set the frame, emission and optionally the emission texture.

```cpp
auto scene = new trace_scene{};             // create a scene
auto environment = add_environment(scene);  // create an environment
environment->frame = identity3x4f;          // set identity transform
auto tex = add_scene(scene, "sky");         // add hdr texture
environment->emission = {1,1,1};            // set emission scale
environment->emission_tex = tex;            // add emission texture
```

## Scene tesselation

The evaluation functions defined above and the ray intersection functions do
not support subdivision surfaces or displaced shapes directly. Instead,
shapes should be converted to indexed meshes using `tesselate_shape(shape)`
for a specific shape, or `tesselate_shapes(scene, progress)` for the
whole scene. Note that tesselations are destructive, meaning that the original
shape data is lost. This is done to avoid copying whenever possible.

```cpp
auto scene = new trace_scene{...};          // create a complete scene
void tesselate_shapes(scene);               // tesselate shapes in the scene
```

## Evaluation of scene properties

Yocto/Trace defines several function to evaluate scene properties.
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
auto texture = scene->textures.front();        // get first texture
auto col = eval_texture(texture,{0.5,0.5});    // eval texture
```

Use `eval_material(material, texcoord)` to evaluate material textures and
combine them with parameter values. The function returns a
`scene_material_sample` that has the same parameters of a material but no
textures defined.

```cpp
auto scene = new trace_scene{...};             // create a complete scene
auto material = scene->materials.front();        // get first material
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
auto instance = scene->instances.front();      // get first instance
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
auto environment = scene->environments.front();  // get first environment
auto envi = eval_environment(environment, dir);  // eval environment
```

## Scene tesselation

The evaluation functions defined above and the ray intersection functions do
not support subdivision surfaces or displaced shapes directly. Instead,
shapes should be converted to indexed meshes using `tesselate_shape(shape)`
for a specific shape, or `tesselate_shapes(scene, progress)` for the
whole scene. Note that tesselations are destructive, meaning that the original
shape data is lost. This is done to avoid copying whenever possible.

```cpp
auto scene = new sceneio_scene{...};          // create a complete scene
void tesselate_shapes(scene);               // tesselate shapes in the scene
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
auto texture = scene->textures.front();        // get first texture
auto col = eval_texture(texture,{0.5,0.5});    // eval texture
```

Use `eval_material(material, texcoord)` to evaluate material textures and
combine them with parameter values. The function returns a
`scene_material_sample` that has the same parameters of a material but no
textures defined.

```cpp
auto scene = new sceneio_scene{...};             // create a complete scene
auto material = scene->materials.front();        // get first material
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
auto instance = scene->instances.front();      // get first instance
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
auto environment = scene->environments.front();  // get first environment
auto envi = eval_environment(environment, dir);  // eval environment
```
