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
with `tesselate_shapes(scene, params, progress)`, then initialize the scene
bvh and lights, with `init_bvh(scene, params, progress)` and
`init_lights(scene, params, progress)`, and then call
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
init_bvh(scene, params, progress);        // init bvh
init_lights(scene, params, progress);     // init lights
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

## Experimental async rendering

The render can run in asynchronous mode where the rendering process is
lunched and runs until completion. This mode is useful when for interactive
viewing or in modeling-while-rendering applications. The API is very minimal
and only controls the rendering process. Use `trace_start(...)` to start
the async renderer and `trace_stop(...)` to stop it. The renderer starts
by rendering a low resolution preview and then proceeds progressively.

The async renderer takes a `trace_state` struct that tracks the rendering
process and contains all data needed by the async renderer. Rendering progress
ifs given by three callbacks. The first two are the rendering callbacks
defined for offline rendeirng, that return progress report and an image buffer
for each sample. In async mode, a further callback is called after each
pixel is rendered.

During rendering, no scenes changes are allowed, and changes to bvh and lights
are not tracked. This is on purpose since it allows for a simple API while
retaining maximum speed. To changes the scene, first stop the render, than
apply scene changes, than update lights and bvh if needed, and finally restart
the renderer.
To update the BVH, use `update_bvh(scene, instances, shapes, params)` where
`instances` and `shapes` are the list of modified ids.

```cpp
auto scene = new trace_scene{...};        // initialize scene
auto params = trace_params{};             // default params
init_bvh(scene, params, progress);        // init bvh
init_lights(scene, params, progress);     // init lights
auto imgprogress = [](const image<vec4f>& render, // image progress
      int sample, int samples) {
  display_image(render);                  // display image
};
auto pxlprogress = [](const image<vec4f>& render, // pixel progress
      int sample, int samples, const vec2i& ij) {
  display_pixel(render, ij);              // display pixel
};
auto state = new trace_state{};           // allocate state
trace_start(state, scene, params, {},     // start async renderer
  imgprogress, pxlprogress);              // function returns immediately
run_app(...);                             // application runs normally here
trace_stop(state);                        // stop async renderer
modify_scene(...);                        // make scene changes
trace_start(state, scene, params, {},     // re-start async renderer
  imgprogress, pxlprogress);              // function returns immediately
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
where `<object>` is the object type name. In these functions,
the name is optional and, if left blank, a unique name will be generated
automatically. For each object type, properties can be set directly.
As a convenience, Yocto/Trace defines several functions to set objects
properties.

For cameras, use `set_frame(camera, frame)` to set the local to world frame,
`set_lens(camera, lens, aspect, film, ortho)` to set the camera projection using
photographic lens parameters, and `set_focus(camera, aperture, focus)` to set
the camera aperture and focus distance.

```cpp
auto scene = new trace_scene{};          // create a scene
auto camera = add_camera(scene);         // create a camera named cam
set_frame(camera,identity3x4f);          // set frame to identity
set_lens(camera,0.050,1.5,0.036);     // set as 50mm lens 3:2 aspect on 35mm
set_aperture(camera,0.01,10);            // set 10mm aperture focused at 10m
```

For instances, use `set_frame(instance, frame)` to set the local to world frame,
`set_shape(instance, shape)` and `set_material(instance, material)`
to set the shape and material pointers. Since adding instances of single
shapes is common in simpler scenes, the function
`add_complete_instance(scene, name)` adds an instance with a new shape and
material.

```cpp
auto scene = new trace_scene{};             // create a scene
auto instance = add_instance(scene);        // create an instance named ist
set_frame(instance,identity3x4f);           // set frame to identity
auto shape = add_shape(scene);
set_shape(instance,shape);                  // set shape pointer
auto material = add_material(scene);
set_material(instance,material);            // set material pointer
auto instance1 = add_complete_instance(scene);  // create an instance
print_info(instance1->shape);                   // with a new shape
print_info(instance1->material);                // and  a new material
```

For textures, use `set_texture(texture, img)` to set the texture
to the specified image. The function has overloads for images with
one or three channels and with float or byte channel types.

```cpp
auto scene = new trace_scene{};             // create a scene
auto texture = add_texture(scene);          // create a texture named tex
set_texture(texture,image<vec4f>{...});     // set as a HDR texture
set_texture(texture,image<vec4b>{...});     // set as a LDR texture
```

For materials, Yocto/Trace defines functions to set each material property.
Each functions take as input the parameter value and an optional texture.
Use `set_emission(material, emission, tex)` to set material emission.
Use `set_color(material, color, tex)` to set the surface color,
`set_specular(material, specular, tex)`,
`set_metallic(material, metallic, tex)`,
`set_transmission(material, transmission, tex)`
for specular, metallic and transmission weights,
`set_ior(material, ior)`, `set_roughness(material, roughness, tex)`,
`set_opacity(material, opacity, tex)` for surface ior, roughness and opacity,
`set_scattering(material, scattering,tex)` for volumetric scattering and
`set_thin(material, thin)` for the thin flag.

```cpp
auto scene = new trace_scene{};               // create a scene
auto matte = add_texture(scene);              // create a matte material
set_color(matte, {1,0.5,0.5}, add_texture(scene)); // textured albedo
auto plastic = add_texture(scene); // create a plastic material
set_color(plastic, {0.5,1,0.5});              // constant color
set_specular(plastic, 1);                     // constant specular
set_roughness(plastic, 0.1, add_texture(scene)); // textured roughness
auto metal = add_texture(scene);              // create a metal material
set_color(metal, {0.5,0.5,1});                // constant color
set_specular(metal, 1);                       // constant specular
set_roughness(metal, 0.1);                    // constant roughness
auto tglass = add_texture(scene);             // create a thin glass material
set_color(tglass, {1,1,1});                   // constant color
set_specular(tglass, 1);                      // constant specular
set_transmission(tglass, 1);                  // constant transmission
auto glass = add_texture(scene);              // create a glass material
set_color(glass, {1,1,1});                    // constant color
set_specular(glass, 1);                       // constant specular
set_transmission(glass, 1);                   // constant transmission
set_thin(glass, false);                       // volumetric material
auto subsurf = add_texture(scene);            // create a subsurface material
set_color(subsurf, {1,1,1});                  // constant color
set_specular(subsurf, 1);                     // constant specular
set_transmission(subsurf, 1);                 // constant transmission
set_thin(subsurf, false);                     // volumetric material
set_scattering(subsurf, {0.5,1,0.5});         // volumetric scattering
```

For shapes, Yocto/Trace defines functions to set shape element indices
and vertex properties. Use `set_points(shape, points)`,
`set_lines(shape, lines)`, `set_triangles(shape, triangles)`, and
`set_quads(shape, quads)` to set indexed meshes indices as points, lines,
triangles and quads respectively.
Use `set_positions(shape, positions)`, `set_normals(shape, normals)`,
`set_texcoords(shape, texcoords)`, `set_colors(shape, colors)`,
`set_radius(shape, radius)`, and `set_tangents(shape, tangents)`
to set positions, normals, texture coordinates, colors, radius and
tangent spaces respectively.

```cpp
auto scene = new trace_scene{};             // create a scene
auto shape = add_shape(scene);              // create a shape named shp
set_triangles(shape, vector<vec3i>{...});   // set triangle indices
set_positions(shape, vector<vec3f>{...});   // set positions
set_normals(shape, vector<vec3f>{...});     // set normals
set_texcoords(shape, vector<vec2f>{...});   // set texture coordinates
```

Use `set_fvquads(shape, quadspos, quadsnorm, quadsuv)` to set the shapes as
a facce-varying quad mesh. Use `set_displacement(shape, disp, tex)` to
set the displacement map and scale and `set_subdivision(shape, level, cc)`
to set the subdivision level and whether to use Catmull-Clark or linear
subdivision.

```cpp
auto scene = new trace_scene{};             // create a scene
auto shape = add_shape(scene);              // create a shape named shp
set_fvquads(shape, vector<vec4i>{...},      // set face-varying indices
  {}, vector<vec4i>{...});                  // for positions and textures
set_positions(shape, vector<vec3f>{...});   // set positions
set_texcoords(shape, vector<vec2f>{...});   // set texture coordinates
set_subdivision(shape, 2, true);            // set Catmull-Clark subdivision
set_displacement(shape, 1, tex);            // sete displacement map tex
```

For environments, use `set_frame(environment, frame)` to set the local to
world frame, `set_emission(instance, emission, emission_tex)` to set the
environment emission and emission texture.

```cpp
auto scene = new trace_scene{};             // create a scene
auto environment = add_environment(scene);  // create an environment
set_frame(environment, identity3x4f);       // set identity transform
auto tex = add_scene(scene, "sky");         // add hdr texture
set_emission(environment, {1,1,1}, tex);    // add emission scale and texture
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
`init_bvh(scene, params, progress)`. The function takes a scene as input and
builds a BVH that is stored internally. The BVH build strategy,
and whether or not Embree is used is determined by the `params` settings.
The function takes also an optional progress callback that is called as
the BVH is build to report build progress. After initialization,
if scene shapes and instances are modified, the BVH can be updated with
`update_bvh(...)`.

```cpp
auto scene = new trace_scene{...};               // create a complete scene
auto params = scene_bvh_params{};                // default params
auto progress = [](const string& message,        // progress callback
   int current, int total) {
  print_info(message, current, total);
};
init_bvh(scene, params, progress);               // build bvh
params.type = scene_bvh_type::embree_default;    // set build type as Embree
init_bvh(scene, params);        // build Embree bvh with no progress report
```

Use `intersect_scene_bvh(scene, ray)` to intersect a ray with a scene,
and `intersect_instance_bvh(instance, ray)` to intersect a ray with an
instance. Both functions return a `scene_intersection` that includes
a hit flag, the instance id, the shape element id, the shape element uv
and intersection distance.

```cpp
auto scene = new trace_scene{...};          // create a complete scene
init_bvh(scene, {});                        // build default bvh
auto ray = ray3f{...};                      // ray
auto isec = intersect_scene_bvh(scene, ray);// ray-scene intersection
if (isec.hit) print_info(isec);
auto instance = scene->instances.first();   // get instance
auto iisec = intersect_instance_bvh(instance, ray); // ray-instance int,
if (iisec.hit) print_info(isec);
```
