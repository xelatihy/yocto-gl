# Yocto/Scene: Scene representation

Yocto/Scene define a simple scene representation, and related utilities,
used in the Yocto/GL path tracer and for scene IO.
Yocto/Scene is implemented in `yocto_scene.h` and `yocto_scene.cpp`.

## Scene representation

Scenes are stored in `scene_model` structs and are comprised of arrays of
cameras, instances, shapes, materials, textures and environments.
The various objects are stored as values in arrays named like the object type.
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

Objects are added to the scene by directly adding elements to the corresponding
arrays. References to elements are expressed as indices to the
corresponding arrays. For each element type, properties can be set directly.
Also, all scene objects are values, so you can work freely with them without
concerning yourself with memory management. The mantra we followed here is
that "if you know how to use `std::vector`, you know how to use scenes".

Here is an sketch of how to create a shape instance in a scene.

```cpp
auto scene = scene_model{};         // create a scene
auto shape = scene_shape{};         // create a shape and add it
set_shape_properties(shape, ...);
scene.shapes.push_back(shape);
scene.materials.push_back({});      // create a black material directly
auto instance = scene_instance{};   // create an instance of last added shape
instance.shape = (int)scene.shapes.size()-1;
instance.material = (int)scene.materials.size()-1;
```

Yocto/Scene defines several function to evaluate scene properties.
Use `compute_bounds(scene)` to compute the scene bounding boxes,  
`scene_stats(scene)` to get scene stats and
`scene_validation(scene)` to validate scene objects.

```cpp
auto scene = scene_scene{...};              // create a complete scene
auto bbox = compute_bounds(scene);          // get bounds
auto stats = scene_stats(scene);            // get stats
for(auto stat : stats) print_info(stat);    // print stats
auto errors = validate_stats(scene);        // get validation errors
for(auto error : errors) print_error(error);// print error
```

## Cameras

Cameras, represented by `scene_camera`, are based on a simple lens model.
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

To create cameras, you should set the camera frame, the camera view,
via lens, aspect and film, and optionally the camera aperture and focus.

```cpp
auto camera = scene_camera{};    // create a camera
camera.frame = identity3x4f;     // set frame to identity
camera.lens = 0.050;             // set as 50mm lens
camera.aspect = 1.5;             // set 3:2 aspect ratio
camera.film = 0.036;             // set the film as 35mm
camera.aperture = 0.01;          // set 10mm aperture
camera.focus = 10;               // set the focus at 10m
```

Use `get_camera(scene, name)` to get a camera by name or the default camera
is the name is not given. Use `eval_camera(camera, image_uv, lens_uv)`
to get a camera ray from the normalized image coordinates `image_uv` and
lens coordinates `lens_uv`.

```cpp
auto scene = scene_model{...};                 // create a complete scene
auto& camera = get_camera(scene);              // get default camera
auto ray = eval_camera(camera,{0.5,0.5},{0,0});// get ray though image center
```

## Instances

Instances, represented as `scene_instance`, place shapes in the scene by
defining their coordinate frame, a shape index and a material index.
Through the use of instancing, Yocto/Scene scales well to large environments
without introducing more complex mechanisms.

For instances, you should set the instance frame, shape and material.

```cpp
auto instance = scene_instance{};    // create an instance
instance.frame = identity3x4f;       // set frame to identity
instance.shape = shape_index;        // set shape index
instance.material = material_index;  // set material index
```

Several functions are defined to evaluate the geometric and material
properties of points on shapes and instances, indicated by the shape element id
and, when needed, the shape element barycentric coordinates. The difference
between the shape and instance methods is that the former returns quantities
in object space, while the latter in world space.
Use `eval_position(...)` to evaluate the point position,
`eval_normal(...)` to evaluate the interpolate point normal,
`eval_texcoord(...)` to evaluate the point texture coordinates,
`eval_element_normal(...)` to evaluate the point geometric normal, and
`eval_color(...)` to evaluate the interpolate point color.
Use `eval_material(...)` as a convenience function to evaluate material
properties of instance points.

```cpp
auto eid = 0; auto euv = vec3f{0.5,0.5};       // element id and uvs
auto pos  = eval_position(instance, eid, euv); // eval point position
auto norm = eval_normal(instance, eid, euv);   // eval point normal
auto st   = eval_texcoord(instance, eid, euv); // eval point texture coords
auto col  = eval_color(instance, eid, euv);    // eval point color
auto gn   = eval_element_normal(instance, eid, euv); // eval geometric normal
auto mat  = eval_material(instance, eid, euv); // eval point material
```

## Environments

Environments, represented as `scene_environment`, store the background
illumination as a scene. Environments have a frame, to rotate illumination,
an emission term and an optional emission texture.
The emission texture is an HDR environment map stored in a LatLon
parametrization.

For environments, set the frame, emission and optionally the emission texture.

```cpp
auto& environment = scene_environment{};  // create an environment
environment.frame = identity3x4f;         // set identity transform
environment.emission = {1,1,1};           // set emission scale
environment.emission_tex = texture_index; // add emission texture
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

## Shapes

Shapes, represented by `scene_shape`, are indexed meshes of elements.
Shapes can contain only one type of element, either
points, lines, triangles or quads. Shape elements are parametrized as in
[Yocto/Geometry](yocto_geometry.md).
Vertex properties are defined as separate arrays and include
positions, normals, texture coords, colors, radius and tangent spaces.
Additionally, Yocto/Scene supports face-varying primitives, as `scene_fvshape`,
where each vertex data has its own topology.

Shapes also work as a standalone mesh representation throughout the
library and can be used even without a scene.

For shapes, you should set the shape elements, i.e. point, limes, triangles
or quads, and the vertex properties, i.e. positions, normals, texture
coordinates, colors and radia. Shapes support only one element type.

```cpp
auto shape = scene_shape{};            // create a shape
shape.triangles = vector<vec3i>{...};  // set triangle indices
shape.positions = vector<vec3f>{...};  // set positions
shape.normals = vector<vec3f>{...};    // set normals
shape.texcoords = vector<vec2f>{...};  // set texture coordinates
```

Several functions are defined to evaluate the geometric properties of points
of shapes, indicated by the shape element id and, when needed, the shape element
barycentric coordinates.
Use `eval_position(...)` to evaluate the point position,
`eval_normal(...)` to evaluate the interpolate point normal,
`eval_texcoord(...)` to evaluate the point texture coordinates,
`eval_element_normal(...)` to evaluate the point geometric normal, and
`eval_color(...)` to evaluate the interpolate point color.

```cpp
auto eid = 0; auto euv = vec3f{0.5,0.5};    // element id and uvs
auto pos  = eval_position(shape, eid, euv); // eval point position
auto norm = eval_normal(shape, eid, euv);   // eval point normal
auto st   = eval_texcoord(shape, eid, euv); // eval point texture coords
auto col  = eval_color(shape, eid, euv);    // eval point color
auto gn   = eval_element_normal(shape, eid, euv); // eval geometric normal
```

Shape support random sampling with a uniform distribution using
`sample_shape(...)` and `sample_shape_cdf(shape)`. Sampling works for lines and
triangles in all cases, while for quad it requires that the elements
are rectangular.

```cpp
auto cdf = sample_shape_cdfd(shape);         // compute the shape CDF
auto points = sample_shape(shape, cdf, num); // sample many points
auto point = sample_shape(shape, cdf,        // sample a single point
  rand1f(rng), rand2f(rng));
```

For shapes, we also support the computation of smooth vertex normals with
`compute_normals(shape)` and converting to and from face-varying representations
with `shape_to_fvshape(shape)` and `fvshape_to_shape(fvshape)`.

## Materials

Materials, represented as `scene_material`, are defined by a material type
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

For materials, we need to specify the material type and color at the minimum.
We can further control the appearance by changing surface roughness, index of
refraction and volumetric properties, when appropriate. Here are some examples.

```cpp
auto matte = scene_material{};           // create a matte material
matte.type = scene_material_type::matte;
matte.color = {1,0.5,0.5};               // with base color and
matte.color_tex = texture_id;            // textured albedo
auto glossy =  scene_material{};         // create a glossy material
glossy.type = scene_material_type::glossy;
glossy.color = {0.5,1,0.5};              // with constant color
glossyv.roughness = 0.1;                 // base roughness and a
glossy.roughness_tex = texture_id;       // roughness texture
auto metallic =  scene_material{};       // create a metallic material
glossy.type = scene_material_type::metallic
metal.color = {0.5,0.5,1};               // constant color
metal.roughness = 0.1;                   // constant roughness
auto tglass = scene_material{};          // create a transparent material
tglass.type = scene_material_type::transparent;
tglass.color = {1,1,1};                  // with constant color
auto glass = scene_material{};           // create a refractive material
glass.type = scene_material_type::transparent;
glass.color = {1,0.9,0.9};               // constant color
auto subsurf = scene_material{};         // create a refractive material
subsurf.type = scene_material_type::subsurface;
subsurf.color = {1,1,1};                 // that transmits all light
subsurf.scattering = {0.5,1,0.5};        // and has volumetric scattering
```

Lights are not explicit in Yocto/Scene but are specified by assigning emissive
materials.

```cpp
auto light = scene_material{};   // create a material
light.color = {0,0,0};           // that does not reflect light
light.emission = {10,10,10};     // but emits it instead
```

Use `eval_material(material, texcoord)` to evaluate material textures and
combine them with parameter values. The function returns a
`material_point` that has the same parameters of a material but no
textures defined.

```cpp
auto mat = eval_material(scene,material,{0.5,0.5}) // eval material
```

## Textures

Textures, represented as `scene_texture`, contains either 8-bit LDR or
32-bit float HDR images with four channels. Textures can be encoded in either
a linear color space or as sRGBs, depending on an internal flag. The use of
float versus byte is just a memory saving feature.

For textures, set the size, the color space, and _either_ the hdr or ldr pixels.

```cpp
auto hdr_texture = scene_texture{};  // create a texture
hdr_texture.width = 512;             // set size
hdr_texture.height = 512;
hdr_texture.linear = true;           // set color space and pixels for an HDR
hdr_texture.pixelsf = vector<vec4f>{...};
auto ldr_texture = scene_texture{};  // create a texture
ldr_texture.width = 512;             // set size
ldr_texture.height = 512;
ldr_texture.linear = false;          // set color space and pixels for an LDR
ldr_texture.pixelsb = vector<vec4b>{...};
```

Use `eval_texture(texture, uv)` to evaluate the texture at specific uvs.
Textures evaluation returns a color in linear color space, regardless of
the texture representation.

```cpp
auto col = eval_texture(texture,{0.5,0.5});   // eval texture
```

## Subdivs

Subdivs, represented as `scene_subdiv`, support tesselation and displacement
mapping. Subdivs are represented as facee-varying shapes.
Subdivs specify a level of subdivision and can be subdivide elements
either linearly or using Catmull-Clark subdivision. Subdivs also support
displacement by specifying both a displacement texture and a displacement amount.
Differently from most systems, in Yocto/Scene displacement is specified
in the shape and not the material. Subdivs only support tesselation to shapes,
but do not directly support additional evaluation of properties.
Subdivs specified to the shape index to which they are subdivided into.

In this case, set the quads for positions, normals and texture coordinates.
Also set the subdivision level, and whether to use Catmull-Clark or linear
subdivision. Finally, displacement can also be applied by setting a displacement
scale and texture.

```cpp
auto subdiv = scene_sundiv{};             // create a subdiv
subdiv.quadspos = vector<vec4i>{...};     // set face-varying indices
subdiv.quadstexcoord = vector<vec4i>{...};// for positions and textures
subdiv.positions = vector<vec3f>{...};    // set positions
subdiv.texcoords = vector<vec2f>{...};    // set texture coordinates
subdiv.subdivisions = 2;                  // set subdivision level
subdiv.catmullclark = true;               // set Catmull-Clark subdivision
subdiv.displacement = 1;                  // set displacement scale
subdiv.displacement_tex = texture_id;     // and displacement map
```

Most properties on subdivs cannot be directly evaluated, nor they are
supported directly in scene processing. Instead, subdivs are converted to
indexed shapes using `tesselate_subdiv(subdiv, shape)` for a specific subdiv,
or `tesselate_subdivs(scene)` for the whole scene.

```cpp
tesselate_subdivs(scene);     // tesselate all subdivs in the scene
```

## Face-Varying shapes

We also support standalone face-varying shapes, that are not stored in the scene
(see subdivs above). In this case, set the quads for positions,
normals and texture coordinates.

```cpp
auto shape = scene_fvshape{};               // create a shape
shape.quadspos = vector<vec4i>{...};        // set face-varying indices
shape.quadstexcoord = vector<vec4i>{...};   // for positions and textures
shape.positions = vector<vec3f>{...};       // set positions
shape.texcoords = vector<vec2f>{...};       // set texture coordinates
```

## Example scenes

Yocto/Scene has a function to create a simple Cornell Box scene for testing.
There are plans to increase support for more test scenes in the future.

```cpp
auto scene = new sceneio_scene{...};         // create a complete scene
make_cornellbox(scene);                      // make cornell box
```

## Procedural shapes

Yocto/Scene has convenience function to create various procedural shapes,
both for testing and for use in shape creation. These are wrappers to the
corresponding functions in [Yocto/Shape](yocto_shape.md), where we maintain
a comprehensive list of all procedural shapes supported.

Procedural shapes take as input the desired shape resolution, the shape scale,
the uv scale, and additional parameters specific to that procedural shape.
These functions return a quad mesh, stored as a `scene_shape` struct.
Use `make_rect(...)` for a rectangle in the XY plane,
`make_bulged_rect(...)` for a bulged rectangle,
`make_recty(...)` for a rectangle in the XZ plane,
`make_bulged_recty(...)` for a bulged rectangle in the XZ plane,
`make_box(...)` for a box,
`make_rounded_box(...)` for a rounded box,
`make_floor(...)` for a floor in the XZ plane,
`make_bent_floor(...)` for a bent floor,
`make_sphere(...)` for a sphere obtained from a cube,
`make_uvsphere(...)` for a sphere tessellated along its uvs,
`make_capped_uvsphere(...)` for a sphere with flipped caps,
`make_disk(...)` for a disk obtained from a quad,
`make_bulged_disk(...)` for a bulged disk,
`make_uvdisk(...)` for a disk tessellated along its uvs,
`make_uvcylinder(...)` for a cylinder tessellated along its uvs,
`make_rounded_uvcylinder(...)` for a rounded cylinder.

```cpp
// make shapes with 32 steps in resolution and scale of 1
auto shape_01 = make_rect({32,32}, {1,1});
auto shape_02 = make_bulged_rect({32,32}, {1,1});
auto shape_03 = make_recty({32,32}, {1,1});
auto shape_04 = make_box({32,32,32}, {1,1,1});
auto shape_05 = make_rounded_box({32,32,32}, {1,1,1});
auto shape_06 = make_floor({32,32}, {10,10});
auto shape_07 = make_bent_floor({32,32}, {10,10});
auto shape_08 = make_sphere(32, 1);
auto shape_09 = make_uvsphere({32,32}, 1);
auto shape_10 = make_capped_uvsphere({32,32}, 1);
auto shape_11 = make_disk(32, 1);
auto shape_12 = make_bulged_disk(32, 1);
auto shape_13 = make_uvdiskm({32,32}, 1);
auto shape_14 = make_uvcylinder({32,32,32}, {1,1});
auto shape_15 = make_rounded_uvcylinder({32,32,32}, {1,1});
```

Yocto/Shape defines a few procedural face-varying shapes with similar interfaces
to the above functions. In this case, the functions return face-varying quads
packed in a `scene_fvshape` struct.
Use `make_fvrect(...)` for a rectangle in the XY plane,
`make_fvbox(...)` for a box,
`make_fvsphere(...)` for a sphere obtained from a cube.

```cpp
// make face-varying shapes with 32 steps in resolution and scale of 1
auto fvshape_01 = make_fvrect({32,32}, {1,1});
auto fvshape_02 = make_fvbox({32,32,32}, {1,1,1});
auto fvshape_03 = make_fvsphere(32, 1);
```

Yocto/Shape provides functions to create predefined shapes helpful in testing.
These functions take only a scale and often provide only the positions as
vertex data. These functions return either triangles, quads, or
face-varying quads in a `scene_shape` or `scene_fvshape` struct.
Use `make_monkey(...)` for the Blender monkey as quads and positions only,
`make_quad(...)` for a simple quad,
`make_quady(...)` for a simple quad in the XZ plane,
`make_cube(...)` for a simple cube as quads and positions only,
`make_fvcube(...)` for a simple face-varying unit cube,
`make_geosphere(...)` for a geodesic sphere as triangles and positions only.
These functions return a `scene_shape` or `scene_fvshape`.

```cpp
auto monkey = make_monkey(1);
auto quad   = make_quad(1);
auto quady  = make_quady(1);
auto cube   = make_cube(1);
auto geosph = make_geosphere(1);
auto fvcube = make_fvcube(1);
```

Yocto/Shape supports the generation of points and lines sets.
Use `make_lines(...)` to create a line set in the XY plane,
`make_points(...)` for a collection of points at the origin,
adn `make_random_points(...)` for a point set randomly placed in a box.
These functions return points or lines, packed in a `scene_shape` struct.

```cpp
auto lines_01 = make_lines({4, 65536},      // line steps and number of lines
                           {1, 1}, {1, 1},  // line set scale and uvscale
                           {0.001, 0.001}); // radius at the bottom and top
// procedural points return points, positions, normals, texcoords, radia
auto [points, positions, normals, texcoords, radius] = make_points(65536);
auto points_01 = make_points(65536,        // number of points
                             1,            // uvscale
                             0.001);       // point radius
auto points_02 = make_random_points(65536, // number of points
                             {1, 1, 1}, 1, // line set scale and uvscale
                             0.001);       // point radius
```

Yocto/Shape also defines a simple functions to generate randomized hairs
on a triangle or quad mesh. Use `make_hair(...)` to create a hair shape
from a triangle and quad mesh, and return a line set.

```cpp
// Make a hair ball around a shape
auto lines =  make_hair(
  make_sphere(),  // sampled surface
  {8, 65536},     // steps: line steps and number of lines
  {0.1, 0.1},     // length: minimum and maximum length
  {0.001, 0.001}, // radius: minimum and maximum radius from base to tip
  {0, 10},        // noise: noise added to hair (strength/scale)
  {0, 128},       // clump: clump added to hair (strength/number)
  {0, 0});        // rotation: rotation added to hair (angle/strength)
```

Finally, Yocto/Shape defines a function to create a quad mesh from a heighfield.
Use `make_heightfield(...)` to create a heightfield meshes.

```cpp
auto heightfield = vctor<float>{...};             // heightfield data
auto shape = make_heightfield(size, heightfield); // make heightfield mesh
```
