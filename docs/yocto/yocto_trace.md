# Yocto/Trace: Path tracing

Yocto/Trace is a simple path tracer written on the Yocto/Scene model.
Yocto/Trace is implemented in `yocto_trace.h` and `yocto_trace.cpp`.

## Rendering features

Yocto/Trace uses the [Yocto/Scene](yocto_scene.md) model to represent
scene data. Consult the scene documentation for a detailed description
of the scene data model. Yocto/Trace supports rendering of shapes
represented as indexed meshes of points, lines, triangles and quads.
The renderer use instancing throughout to gain scalability with large
environment, while still remaining easy to use. Materials are modeled
similarly to the Disney's BSDF with integrated subsurface scattering.
Lights are not specified explicitly but derived by the renderer
from emissive surface and environment maps in the scene.

Yocto/Trace implements a fast but simple path tracer. The algorithm
uses relatively advanced sampling features to reduce the numbers of
rays used to render an image. At the same time, we focus on general
solutions that work well in all cases, instead of tuning the algorithm
to work well in some specific cases.

Yocto/Trace is a progressively path tracer that computes images
a few samples at a time. This is true for both offline and online
rendering. This design allows for partial images to be used to
give feedback before the render is finished. In the future, this
might allow for more rendering features like adaptive rendering
or stop-and-resume renders.

## Rendering a scene

To render a scene, first tesselate shapes for subdivs and displacement,
with `tesselate_shapes(scene, params, progress)`, then call
`trace_image(scene, camera, params, progress, image_progress)`.
In these functions, `params` are the rendering options used to
render the scene, `progress` is a callback function used to
report rendering progress and `image_progress` is a callback
function used to report partial images during rendering.
Both callback functions are optional.

```cpp
auto scene = new trace_scene{...};        // initialize scene
auto params = trace_params{};             // default params
auto progress = [](const string& message, // progress callback
      int current, int total) {
  print_info(message, current, total);
};
tesselate_shapes(scene, params, progress);// tesselate shapes if needed
trace_image(scene, params, progress);     // render image
auto improgress = [](const image<vec4f>& render, // image progress
      int sample, int samples) {
  if(sample % 16 == 0)                    // every 16 samples
    save_image(filename, render);         // save image
};
trace_image(scene, params, progress, improgress); // render image
```

## Rendering options

The rendering process is configured with the `trace_params` settings.
`sampler` determines the algorithm used for rendering: `path` is the default
algorithm and probably the one you want to use, `naive` is a simpler path
tracing that may be used for testing, `eyelight` produces quick previews
of the screen geometry, `falsecolor` is a debug feature to view scenes
according to the `falsecolor` setting, and `albedo` and `normal` are
used to produce denoising buffers.

THe image resolution is set by `resolution` and measures the resolution
of the longest axis. `samples` is the number of per-pixel samples
used while rendering and is the only parameter used to control the
tradeoff between noise and speed. `bounces` is the maximum number of bounces
and should be high for scenes with glass and volumes, but otherwise a low
number would suffice.

The remaining parameters are approximation used to reduce noise, at the
expenses of bias. `clamp` remove high-energy fireflies. `nocaustics` removes
certain path that cause caustics. `tentfilter` apply a linear filter to the
image pixels. `envhidden` removes the environment map from the camera rays.

Finally, the `bvh` parameter controls the heuristic used to build the Bvh
and whether the Bvh uses Embree. Please see the description in Yocto/Scene.

`trace_sampler_names`, `trace_falsecolor_names` and `trace_bvh_names`
define string names for various enum values that can used for UIs or CLIs.

```cpp
// high quality rendering
auto hq_params = trace_params{};             // default params
hq_params.sampler    = trace_sampler_type::path; // path tracing
hq_params.resolution = 1280;                 // high-res render
hq_params.samples    = 1024;                 // high-sample count
hq_params.bounce     = 64;                   // high max bounces
// geometry previewing
auto pp_params = trace_params{};             // default params
pp_params.sampler    = trace_sampler_type::eyelight; // geometry preview
pp_params.resolution = 720;                  // medium-res render
pp_params.samples    = 16;                   // low-sample count
// scene debug viewing
auto db_params = trace_params{};             // default params
db_params.sampler    = trace_sampler_type::falsecolor; // debug false colors
db_params.falsecolor = trace_falsecolor_type::normal; // normal viewing
db_params.resolution = 720;                  // medium-res render
db_params.samples    = 16;                   // low-sample count
```

## Low-Level Rendering API

Yocto/Trace supports a low-level rendering API that is more flexible for
applications that need it. In this modality, you have to manage all state
manually.

Yocto/Trace internally uses a `trace_bvh` BVH for ray-scene intersection,
a `trace_lights` objects to store lighting information, and a `trace_state`
object to tracks the rendering process and contains all data needed to perform
progressive computation. Each object need to be initialized separately.

To render a scene, first tesselate shapes for subdivs and displacement,
with `tesselate_shapes(scene, params, progress)`, then initialize the scene
bvh and lights, with `init_bvh(bvh, scene, params, progress)` and
`init_lights(lights, scene, params, progress)`, then initialize the state
with `init_state(state, scene, params)` and then call
`trace_image(state, scene, camera, lights, params, progress, image_progress)`.
In these functions, `params` are the rendering options used to
render the scene, `progress` is a callback function used to
report rendering progress and `image_progress` is a callback
function used to report partial images during rendering.
Both callback functions are optional.

```cpp
auto scene = new trace_scene{...};        // initialize scene
auto params = trace_params{};             // default params
auto progress = [](const string& message, // progress callback
      int current, int total) {
  print_info(message, current, total);
};
tesselate_shapes(scene, params, progress);// tesselate shapes if needed
auto bvh = new trace_bvh{};                   // trace bvh
init_bvh(bvh, scene, params, progress);       // init bvh
auto lights = new trace_lights{};             // trace lights
init_lights(lights, scene, params, progress); // init lights
auto state = new trace_state{};               // trace state
init_state(state, scene, params, progress);   // init lights
trace_image(state, scene, camera, bvh,       // render image
  lights, params, progress);
auto improgress = [](const image<vec4f>& render, // image progress
      int sample, int samples) {
  if(sample % 16 == 0)                        // every 16 samples
    save_image(filename, render);             // save image
};
trace_image(state, scene, camera, hvh,        // render image
  lights, params, progress, improgress);
```

## Experimental async rendering

The render can run in asynchronous mode where the rendering process is
lunched and runs until completion. This mode is useful when for interactive
viewing or in modeling-while-rendering applications. The API is very minimal
and only controls the rendering process. Use `trace_start(...)` to start
the async renderer and `trace_stop(...)` to stop it. The renderer starts
by rendering a low resolution preview and then proceeds progressively.

The async renderer takes a `trace_state` struct that tracks the rendering
process and contains all data needed by the async renderer. Rendering progress
is given by three callbacks. The first two are the rendering callbacks
defined for offline rendeirng, that return progress report and an image buffer
for each sample. In async mode, a further callback is called after each
pixel is rendered.

During rendering, no scenes changes are allowed, and changes to bvh and lights
are not tracked. This is on purpose since it allows for a simple API while
retaining maximum speed. To change the scene, first stop the render, than
apply scene changes, than update lights and bvh if needed, and finally restart
the renderer.
To update the BVH, use `update_bvh(scene, instances, shapes, params)` where
`instances` and `shapes` are the list of modified ids.

```cpp
auto scene = new trace_scene{...};            // initialize scene
auto params = trace_params{};                 // default params
auto bvh = new trace_bvh{};                   // trace bvh
init_bvh(bvh, scene, params, progress);       // init bvh
auto lights = new trace_lights{};             // trace lights
init_lights(light, scene, params, progress);  // init lights
auto imgprogress = [](const image<vec4f>& render, // image progress
      int sample, int samples) {
  display_image(render);                      // display image
};
auto pxlprogress = [](const image<vec4f>& render, // pixel progress
      int sample, int samples, const vec2i& ij) {
  display_pixel(render, ij);                  // display pixel
};
auto state = new trace_state{};               // allocate state
trace_start(state, scene, camera, bvh,        // start async renderer
  lights, params, {},                         // and return immediately
  imgprogress, pxlprogress);                  // communicates via callbacks
run_app(...);                                 // application continues
trace_stop(state);                            // stop async renderer
modify_scene(...);                            // make scene changes
trace_start(state, scene, camera, bvh,        // re-start async renderer
  lights, params, {},                         // return immediately
  imgprogress, pxlprogress);                  // communicates via callbacks
```

## Scene representation

Scenes are stored in `trace_scene` structs and are comprised of array
of objects whose memory is owned by the scene.
Scenes are comprised of camera, instances, shapes, materials, textures
and environments. The scene representation is geared toward modeling
physically-based environments and is similar to the
[Yocto/Trace](yocto_sceneio.md) scene.
In Yocto/Trace, lights are not explicitly defined, but
implicitly comprised of instances with emissive materials and environment maps.
All scenes and objects properties are accessible directly.
Some scene data, like lights and ray-acceleration structures,
are computed when necessary and should not be accessed directly.

Cameras, instances and environments have coordinate frames to define
the local to world transformation.
Frames are presented as affine 3x4 matrices and are intended to be
rigid transforms, although most scene processing support frames with
scaling.

**Cameras**, represented by `trace_camera`, are based on a simple lens model.
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

**Textures**, represented as `trace_texture` contain either 8-bit LDR or
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
`trace_shape` type. Shapes can contain only one type of element, either
points, lines, triangles or quads. Shape elements are parametrized as in
[Yocto/Geometry](yocto_geometry.md).
Vertex properties are defined as separate arrays and include
positions, normals, texture coords, colors, radius and tangent spaces.
Additionally, Yocto/Trace supports face-varying primitives where
each vertex data has its own topology.

Shapes supports tesselation and displacement mapping. Shape specify a
level of subdivision and can be subdivide elements either linearly
or using Catmull-Clark subdivision. Shapes also support displacement
by specifying both a displacement texture and a displacement amount.
Differently from most systems, in Yocto/Trace displacement is specified
in the shape and not the material, that only controls shading.
Subdivision and displacement are only specified in shapes, but not
taken into account when evaluating shape properties. For this to happen,
shapes have to be tessellated, as shown later.

Shapes are placed in the scene by defining shape **instances** that
take a coordinate frame, a shape pointer and a material pointer.
Instances are represented by the `trace_instance` type.
Instances are represented as `trace_instance` objects. Thought the
use of instancing Yocto/Trace scales well to large environments without
introducing more complex mechanisms.

Scenes might be lit by background illumination defined by **environments**,
represented by the `trace_environment` type. Environments have a frame,
to rotate illumination, an emission term and an optional emission texture.
The emission texture is an HDR environment map stored in a LatLon
parametrization.

## Scene Creation

Objects are added to the scene via `add_<object>(scene)` functions,
where `<object>` is the object type name. For each object type,
properties can be set directly.

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

## Ray-scene intersection

Yocto/Trace supports ray-scene intersection queries accelerated by a two-level
BVH. We provide both our implementation and an Embree wrapper. To perform
ray intersection queries, first initialize the BVH with
`init_bvh(bvh, scene, params, progress)`. The function takes a scene as input
and builds a BVH that is stored internally. The BVH build strategy,
and whether or not Embree is used is determined by the `params` settings.
The function takes also an optional progress callback that is called as
the BVH is build to report build progress. After initialization,
if scene shapes and instances are modified, the BVH can be updated with
`update_bvh(...)`.

```cpp
auto scene = new trace_scene{...};               // create a complete scene
auto params = trace_params{};                    // default params
auto progress = [](const string& message,        // progress callback
   int current, int total) {
  print_info(message, current, total);
};
init_bvh(bvh, scene, params, progress);          // build bvh
params.type = trace_bvh_type::embree_default;    // set build type as Embree
init_bvh(bvh, scene, params);        // build Embree bvh with no progress report
```

Use `intersect_bvh(bvh, ray)` to intersect a ray with a scene,
and `intersect_bvh(bvh, instance, ray)` to intersect a ray with an
instance. Both functions return a `bvh_intersection` that includes
a hit flag, the instance id, the shape element id, the shape element uv
and intersection distance.

```cpp
auto scene = new trace_scene{...};          // create a complete scene
auto bvh = new trace_bvh{};                 // trace bvh
init_bvh(bvh, scene, {});                   // build default bvh
auto ray = ray3f{...};                      // ray
auto isec = intersect_bvh(bvh, ray);        // ray-scene intersection
if (isec.hit) print_info(isec);
auto iisec = intersect_bvh(scene, 0, ray);  // ray-instance intersection
if (iisec.hit) print_info(isec);
```
