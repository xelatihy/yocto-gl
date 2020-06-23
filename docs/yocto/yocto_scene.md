# Yocto/Scene: Scene representation

Yocto/Scene define a simple scene representation, and related utilities,
used to write other libraries in Yocto/GL.
Yocto/Scene is implemented in `yocto_scene.h` and `yocto_scene.cpp`.

## Scene representation

Scenes are stored in `scene_model` structs and are comprised of array
of objects whose memory is owned by the scene.
Scenes are comprised of camera, instances, shapes, materials, textures
and environments. Animation is not currently supported.
The scene representation is geared toward modeling physically-based environments. In Yocto/Scene, lights are not explicitly defined, but
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

_Cameras_, represented by `scene_camera`, are based on a simple lens model.
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

_Textures_, represented as `scene_texture` contain either 8-bit LDR or
32-bit float HDR images of either scalar or color (RGB) type.
HDR images are encoded in linear color space, while LDRs are encoded in sRGB.

_Materials_

// Material for surfaces, lines and triangles.
// For surfaces, uses a microfacet model with thin sheet transmission.
// The model is based on OBJ, but contains glTF compatibility.
// For the documentation on the values, please see the OBJ format.
struct scene_material {
// material data
string name = "";

// material
vec3f emission = {0, 0, 0};
vec3f color = {0, 0, 0};
float specular = 0;
float roughness = 0;
float metallic = 0;
float ior = 1.5;
vec3f spectint = {1, 1, 1};
float coat = 0;
float transmission = 0;
float translucency = 0;
vec3f scattering = {0, 0, 0};
float scanisotropy = 0;
float trdepth = 0.01;
float opacity = 1;
bool thin = true;

// textures
scene_texture* emission_tex = nullptr;
scene_texture* color_tex = nullptr;
scene_texture* specular_tex = nullptr;
scene_texture* metallic_tex = nullptr;
scene_texture* roughness_tex = nullptr;
scene_texture* transmission_tex = nullptr;
scene_texture* translucency_tex = nullptr;
scene_texture* spectint_tex = nullptr;
scene_texture* scattering_tex = nullptr;
scene_texture* coat_tex = nullptr;
scene_texture* opacity_tex = nullptr;
scene_texture* normal_tex = nullptr;
};

_Shapes_ are represented as indexed meshes of elements using the
`scene_shape` type. Shapes can contain only one type of element, either
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
Differently from most systems, in Yocto/Scene displacement is specified
in the shape and not the material, that only controls shading.
Subdivision and displacement are only specified in shapes, but not
taken into account when evaluating shape properties. For this to happen,
shapes have to be tessellated, as shown later.

Shapes are placed in the scene by defining shape _instances_ that
take a coordinate frame, a shape pointer and a material pointer.
Instances are represented by the `scene_instance` type.
Instances are represented as `scene_instance` objects. Thought the
use of instancing Yocto/Scene scales well to large environments without
introducing more complex mechanisms.

Scenes might be lit by background illumination defined by _environments_,
represented by the `scene_environment` type. Environments have a frame,
to rotate illumination, an emission term and an optional emission texture.
The emission texture is an HDR environment map stored in a LatLon
parametrization.

## Scene Creation

Objects are added to the scene via `add_XXX(scene,name)` functions,
where `XXX` is the object type name. In these functions,
the name is optional and, if left blank, a unique name will be generated
automatically. For each object type, properties can be set directly.
As a convenience, Yocto/Scene defines several functions to set objects
properties.

For cameras, use `set_frame(camera,frame)` to set the local to world frame,
`set_lens(camera,lens,aspect,film,ortho)` to set the camera projection using
photographic lens parameters, and `set_focus(camera,aperture,focus)` to set
the camera aperture and focus distance.

```cpp
auto scene = new scene_model{};          // create a scene
auto camera = add_camera(scene, "cam");  // create a camera named cam
set_frame(camera,identity3x4f);          // set frame to identity
set_lens(camera,0.050,3.0/2.0,0.036);    // set as 50mm lens 3:2 aspect on 35mm
set_aperture(camera,0.01,10);            // set 10mm aperture focused at 10m
```

For instances, use `set_frame(instance,frame)` to set the local to world frame,
`set_shape(instance,shape)` and `set_material(instance,material)`
to set the shape and material pointers. Since adding instances of single
shapes is common in simpler scenes, the function
`add_complete_instance(scene,name)` adds an instance with a new shape and
material.

```cpp
auto scene = new scene_model{};             // create a scene
auto instance = add_instance(scene, "ist"); // create an instance named ist
set_frame(instance,identity3x4f);           // set frame to identity
auto shape = add_shape(scene, "shp");
set_shape(instance,shape);                  // set shape pointer
auto material = add_material(scene, "mat");
set_material(instance,material);            // set material pointer

auto instance1 = add_complete_instance(scene, "obj");  // create an instance
print_info(instance1->shape);                          // with a new shape
print_info(instance1->material);                       // and  a new material
```

For textures, use `set_texture(texture, img)` to set the texture
to the specified image. The function has overloads for images with
one or three channels and with float or byte channel types.

```cpp
auto scene = new scene_model{};             // create a scene
auto texture = add_texture(scene, "tex");   // create a texture named tex
set_texture(texture,image<vec3f>{...});     // set as a color HDR texture
set_texture(texture,image<vec3b>{...});     // set as a color LDR texture
set_texture(texture,image<float>{...});     // set as a scalar HDR texture
set_texture(texture,image<byte >{...});     // set as a scalar LDR texture
```

// material properties
void set_emission(scene_material* material, const vec3f& emission,
scene_texture* emission_tex = nullptr);
void set_color(scene_material* material, const vec3f& color,
scene_texture* color_tex = nullptr);
void set_specular(scene_material* material, float specular = 1,
scene_texture* specular_tex = nullptr);
void set_ior(scene_material* material, float ior);
void set_metallic(scene_material* material, float metallic,
scene_texture* metallic_tex = nullptr);
void set_transmission(scene_material* material, float transmission, bool thin,
float trdepth, scene_texture* transmission_tex = nullptr);
void set_translucency(scene_material* material, float translucency, bool thin,
float trdepth, scene_texture* translucency_tex = nullptr);
void set_roughness(scene_material* material, float roughness,
scene_texture* roughness_tex = nullptr);
void set_opacity(scene_material* material, float opacity,
scene_texture* opacity_tex = nullptr);
void set_thin(scene_material* material, bool thin);
void set_scattering(scene_material* material, const vec3f& scattering,
float scanisotropy, scene_texture* scattering_tex = nullptr);
void set_normalmap(scene_material* material, scene_texture* normal_tex);

For shapes, Yocto/Scene defines functions to set shape element indices and vertex properties. Use `set_points(shape, points)`,
`set_lines(shape, lines)`, `set_triangles(shape, triangles)`, and
`set_quads(shape, quads)` to set indexed meshes indices as points, lines,
triangles and quads respectively.
Use `set_positions(shape, positions)`, `set_normals(shape, normals)`,`set_texcoords(shape, texcoords)`, `set_colors(shape, colors)`,
`set_radius(shape, radius)`, and `set_tangents(shape, tangents)`
to set positions, normals, texture coordinates, colors, radius and
tangent spaces respectively.

```cpp
auto scene = new scene_model{};             // create a scene
auto shape = add_shape(scene, "shp");       // create a shape named shp
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
auto scene = new scene_model{};             // create a scene
auto shape = add_shape(scene, "shp");       // create a shape named shp
set_fvquads(shape, vector<vec4i>{...},      // set face-varying indices
  {}, vector<vec4i>{...});                  // for positions and textures
set_positions(shape, vector<vec3f>{...});   // set positions
set_texcoords(shape, vector<vec2f>{...});   // set texture coordinates
set_subdivision(shape, 2, true);            // set Catmull-Clark subdivision
set_displacement(shape, 1, tex);            // sete displacement map tex
```

For environments, use `set_frame(environment, frame)` to set the local to
world frame, `set_emission(instance,emission,emission_tex)` to set the
environment emission and emission texture.

```cpp
auto scene = new scene_model{};             // create a scene
auto environment = add_environment(scene, "env"); // create an environment
set_frame(environment, identity3x4f);       // set identity transform
auto tex = add_scene(scene, "sky");         // add hdr texture
set_emission(environment, {1,1,1}, tex);    // add emission scale and texture
```

// add missing elements
void add_cameras(scene_model* scene);
void add_radius(scene_model* scene, float radius = 0.001f);
void add_materials(scene_model* scene);
void add_sky(scene_model* scene, float sun_angle = pif / 4);

// Trim all unused memory
void trim_memory(scene_model\* scene);

// Clone a scene
void clone_scene(scene_model* dest, const scene_model* scene);

// compute scene bounds
bbox3f compute_bounds(const scene_model\* scene);

// get named camera or default if name is empty
scene_camera* get_camera(const scene_model* scene, const string& name = "");

## Evaluation of scene properties

// Generates a ray from a camera.
ray3f eval_camera(
const scene_camera\* camera, const vec2f& image_uv, const vec2f& lens_uv);

// Evaluates a texture
vec2i texture_size(const scene_texture* texture);
vec3f lookup_texture(
const scene_texture* texture, const vec2i& ij, bool ldr_as_linear = false);
vec3f eval_texture(const scene_texture\* texture, const vec2f& uv,
bool ldr_as_linear = false, bool no_interpolation = false,
bool clamp_to_edge = false);

// Evaluate instance properties
vec3f eval_position(
const scene_instance* instance, int element, const vec2f& uv);
vec3f eval_element_normal(const scene_instance* instance, int element);
vec3f eval_normal(const scene_instance* instance, int element, const vec2f& uv);
vec2f eval_texcoord(
const scene_instance* instance, int element, const vec2f& uv);
pair<vec3f, vec3f> eval_element_tangents(
const scene_instance* instance, int element);
vec3f eval_normalmap(
const scene_instance* instance, int element, const vec2f& uv);
vec3f eval_shading_normal(const scene_instance* instance, int element,
const vec2f& uv, const vec3f& outgoing);
vec3f eval_color(const scene_instance* instance, int element, const vec2f& uv);

// Environment
vec3f eval_environment(
const scene_environment* environment, const vec3f& direction);
vec3f eval_environment(const scene_model* scene, const vec3f& direction);

// Material sample
struct scene_material_sample {
vec3f emission = {0, 0, 0};
vec3f color = {0, 0, 0};
float specular = 0;
float roughness = 0;
float metallic = 0;
float ior = 1.5;
vec3f spectint = {1, 1, 1};
float coat = 0;
float transmission = 0;
float translucency = 0;
vec3f scattering = {0, 0, 0};
float scanisotropy = 0;
float trdepth = 0.01;
float opacity = 1;
bool thin = true;
vec3f normalmap = {0, 0, 1};
};

// Evaluates material and textures
scene_material_sample eval_material(
const scene_material\* material, const vec2f& texcoord);

// Material Bsdf parameters
struct scene_bsdf {
// brdf lobes
vec3f diffuse = {0, 0, 0};
vec3f specular = {0, 0, 0};
vec3f metal = {0, 0, 0};
vec3f coat = {0, 0, 0};
vec3f transmission = {0, 0, 0};
vec3f translucency = {0, 0, 0};
vec3f refraction = {0, 0, 0};
float roughness = 0;
float ior = 1;
vec3f meta = {0, 0, 0};
vec3f metak = {0, 0, 0};
// weights
float diffuse_pdf = 0;
float specular_pdf = 0;
float metal_pdf = 0;
float coat_pdf = 0;
float transmission_pdf = 0;
float translucency_pdf = 0;
float refraction_pdf = 0;
};

// Eval material to obtain emission, brdf and opacity.
vec3f eval_emission(const scene_instance* instance, int element,
const vec2f& uv, const vec3f& normal, const vec3f& outgoing);
// Eval material to obatain emission, brdf and opacity.
scene_bsdf eval_bsdf(const scene_instance* instance, int element,
const vec2f& uv, const vec3f& normal, const vec3f& outgoing);
float eval_opacity(const scene_instance\* instance, int element, const vec2f& uv,
const vec3f& normal, const vec3f& outgoing);
// check if a brdf is a delta
bool is_delta(const scene_bsdf& bsdf);

// Material volume parameters
struct scene_vsdf {
vec3f density = {0, 0, 0};
vec3f scatter = {0, 0, 0};
float anisotropy = 0;
};

// check if we have a volume
bool has_volume(const scene_instance* instance);
// evaluate volume
scene_vsdf eval_vsdf(
const scene_instance* instance, int element, const vec2f& uv);

## Ray-sene intersection

// Strategy used to build the bvh
enum struct scene*bvh_type {
default*,
highquality,
middle,
balanced,
#ifdef YOCTO_EMBREE
embree_default,
embree_highquality,
embree_compact // only for copy interface
#endif
};

// Params for scene bvh build
struct scene*bvh_params {
scene_bvh_type bvh = scene_bvh_type::default*;
bool noparallel = false;
};

// Progress callback called when loading.
using progress_callback =
function<void(const string& message, int current, int total)>;

// Build the bvh acceleration structure.
void init_bvh(scene_model\* scene, const scene_bvh_params& params,
progress_callback progress_cb = {});

// Refit bvh data
void update_bvh(scene_model* scene,
const vector<scene_instance*>& updated_objects,
const vector<scene_shape\*>& updated_shapes, const scene_bvh_params& params);

// Results of intersect functions that include hit flag, the instance id,
// the shape element id, the shape element uv and intersection distance.
// Results values are set only if hit is true.
struct scene_intersection {
int instance = -1;
int element = -1;
vec2f uv = {0, 0};
float distance = 0;
bool hit = false;
};

// Intersect ray with a bvh returning either the first or any intersection
// depending on `find_any`. Returns the ray distance , the instance id,
// the shape element index and the element barycentric coordinates.
scene_intersection intersect_scene_bvh(const scene_model* scene,
const ray3f& ray, bool find_any = false, bool non_rigid_frames = true);
scene_intersection intersect_instance_bvh(const scene_instance* instance,
const ray3f& ray, bool find_any = false, bool non_rigid_frames = true);

## Example scenes

// Make Cornell Box scene
void make_cornellbox(scene_model\* scene);

## Scene stats and validation

// Return scene statistics as list of strings.
vector<string> scene_stats(const scene_model* scene, bool verbose = false);
// Return validation errors as list of strings.
vector<string> scene_validation(
const scene_model* scene, bool notextures = false);

## Scene tesselation

// Apply subdivision and displacement rules.
void tesselate_shapes(scene_model* scene, progress_callback progress_cb = {});
void tesselate_shape(scene_shape* shape);

```

```
