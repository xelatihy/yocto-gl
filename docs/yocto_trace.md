# Yocto/Trace

Path tracer implementation for support for textured mesh
lights, GGX/Phong materials, environment mapping. The interface supports
progressive parallel execution with any splitting strategy, by
generating per-sample random number sequences or fully deterministic
hash-based sampling.

The raytraced scene is a list of instances of basic shapes. Each shape
is a collection of points, lines or triangles with associated normals.
Shapes are instanced by creating instances with soecific local-to-world
trasforms. Instancing shares memory so large scenes can be created easily.
Shaope data is in fact shared with the application and not copied
internally.

Materials are represented as sums of an emission term, a diffuse term and
a specular microfacet term (GGX or Phong). Only opaque for now. We pick
a proper material type for each shape element type (points, lines,
triangles).

Lights are defined as any shape with a material emission term. Additionally
one can also add an environment map. But even if you can, you might want to
add a large triangle mesh with inward normals instead. The latter is more
general (you can even more an arbitrary shape sun). For now only the first
env is used.

We generate our own random numbers guarantying that there is one random
sequence per path. This means you can run the path tracer in any order
serially or in parallel.

For now, we support a straightforward path tracer with explicit direct
illumination using MIS. Also added simpler shaders for a quick preview
and a direct-only renderer.

The library can support raytracing either by building an internal
acceleration structure with Yocto/Bvh or with user supplied intersection
routines for custom intersection.

This library depends in yocto_math.h and yocto_utils.h/.cpp for concurrency.
Eventually the concurrency calls will be move to std functions when
more readily available.
Optionally depend on yocto_bvh.h/.cpp for internal acceleration.
Disable this by setting YTRACE_NO_BVH.


## Usage for Scene Creation

1. create a scene with `make_scene()`
2. add cameras with `add_camera()`, `set_camera()`
3. add add texture with `add_texture()`
4. create material with `add_XXX_material()`
5. add shapes with `add_XXX_shape()`
6. add instances with `add_instance()`
7. add environment maps with `add_environment()`

## Usage for Rendering

1. either build the ray-tracing acceleration structure with
  `init_intersection()` or supply your own with
  `set_intersection_callbacks()`
2. if desired, add logging with `set_logging_callbacks()`
3. prepare lights for rendering `init_lights()`
4. define rendering params with the `trace_params` structure
5. render blocks of samples with `trace_block()`

## Usage for Progressive Rendering

1. either build the ray-tracing acceleration structure with
  `init_intersection()` or supply your own with
  `set_intersection_callbacks()`
2. if desired, add logging with `set_logging_callbacks()`
3. prepare lights for rendering `init_lights()`
4. define rendering params with the `trace_params` structure
5. initoialize the prograssive rendering state with `init_state()`
6. either render sames successively with `trace_next_samples()`
   or starts an asynchronousn renderer
7. get the rendered image with `get_traced_image()`


## History

- v 0.26: thin glass material
- v 0.25: added refraction (still buggy in some cases)
- v 0.24: corrected transaprency bug
- v 0.23: simpler logging
- v 0.22: added additional buffers
- v 0.21: added filters
- v 0.20: state-based api
- v 0.19: explicit material models
- v 0.18: simpler texture creation functions
- v 0.17: move to rgba per-vertex color
- v 0.16: use yocto_math in the interface and remove inline compilation
- v 0.15: move to add api
- v 1.16: internally use yocto_bvh if desired
- v 1.15: added gltf/generic material properties (deprecate old interface)
- v 1.14: normal mapping
- v 1.13: simpler Fresnel handling
- v 1.12: significantly better path tracing
- v 1.11: add progressive sampling to rendering params
- v 0.10: switch to .h/.cpp pair
- v 0.9: doxygen comments
- v 0.8: opaque API (allows for changing internals without altering API)
- v 0.7: internally use pointers for performance transparency
- v 0.6: minor API change for blocks
- v 0.5: [major API change] move to modern C++ interface
- v 0.4: C++ API
- v 0.3: removal of C interface
- v 0.2: use of STL containers
- v 0.1: C++ implementation
- v 0.0: initial release in C99

## Namespace ytrace

Path tracer with MIS support.

### Struct scene

~~~ .cpp
struct scene;
~~~

Trace scene.

### Function make_scene()

~~~ .cpp
scene* make_scene();
~~~

Initialize the scene with the proper number of objects.

### Function free_scene()

~~~ .cpp
void free_scene(scene*& scn);
~~~

Free scene.

### Function add_camera()

~~~ .cpp
int add_camera(scene* scn, const ym::frame3f& frame, float yfov, float aspect,
    float aperture = 0, float focus = 1);
~~~

Adds a camera in the scene.

- Parameters:
    - scn: scene
    - frame: local-to-world frame (x, y, z, o in column major order)
    - yfov: field of view
    - aspect: aspect ratio
    - aperture: lens aperture
    - focus: focus plane distance (cannot be zero)
- Returns:
    - camera id

### Function set_camera()

~~~ .cpp
void set_camera(scene* scn, int cid, const ym::frame3f& frame, float yfov,
    float aspect, float aperture = 0, float focus = 1);
~~~

Sets a camera in the scene.

- Parameters:
    - scn: scene
    - cid: camera id
    - frame: local-to-world frame (x, y, z, o in column major order)
    - yfov: field of view
    - aspect: aspect ratio
    - aperture: lens aperture
    - focus: focus plane distance (cannot be zero)

### Function add_texture()

~~~ .cpp
int add_texture(scene* scn, int width, int height, const ym::vec4f* hdr);
~~~

Adds a texture in the scene.

- Parameters:
    - scn: scene
    - tid: texture id
    - width: width
    - height: height
    - hdr: hdr pixels
- Returns:
    - texture id

### Function add_texture()

~~~ .cpp
int add_texture(scene* scn, int width, int height, const ym::vec4b* ldr);
~~~

Sets a texture in the scene.

- Parameters:
    - scn: scene
    - tid: texture id
    - width: width
    - height: height
    - ldr: ldr pixels (sRGB)
- Returns:
    - texture id

### Function add_texture()

~~~ .cpp
inline int add_texture(scene* scn, const ym::image4f* img);
~~~

Adds a texture in the scene.

- Parameters:
    - scn: scene
    - hdr: hdr image
- Returns:
    - texture id

### Function add_texture()

~~~ .cpp
inline int add_texture(scene* scn, const ym::image4b* img);
~~~

Sets a texture in the scene.

- Parameters:
    - scn: scene
    - ldr: ldr image (sRGB)
- Returns:
    - texture id

### Function add_material()

~~~ .cpp
int add_material(scene* scn);
~~~

Adds a black material to the scene. Use set_material_XXX() functions to
customize it.

- Parameters:
    - scn: scene
- Returns:
    - material id

### Function set_material_emission()

~~~ .cpp
void set_material_emission(
    scene* scn, int mid, const ym::vec3f& ke, int ke_txt);
~~~

Sets the material emission.

- Parameters:
    - scn: scene
    - mid: material id
    - ke: emission, term
    - ke_txt: emission texture (-1 for none)

### Function set_material_normal()

~~~ .cpp
void set_material_normal(scene* scn, int mid, int norm_txt, float scale = 1);
~~~

Sets the material normal map.

- Parameters:
    - scn: scene
    - mid: material id
    - norm_txt: normal map (-1 for none)
    - scale: normal scale

### Function set_material_occlusion()

~~~ .cpp
void set_material_occlusion(scene* scn, int mid, int occ_txt, float scale = 1);
~~~

Sets the material normal map.

- Parameters:
    - scn: scene
    - mid: material id
    - occ_txt: occlusion map (-1 for none)
    - scale: occlusion scale

### Function set_material_microfacet()

~~~ .cpp
void set_material_microfacet(scene* scn, int mid, const ym::vec3f& kd,
    const ym::vec3f& ks, const ym::vec3f& kt, float rs, float op, int kd_txt,
    int ks_txt, int kt_txt, int rs_txt, int op_txt, bool use_phong = false);
~~~

Sets a material reflectance as a microfacet model in the scene.

- Parameters:
    - scn: scene
    - mid: material id
    - kd: diffuse term
    - ks: specular term
    - kt: transmission term
    - rs: specular roughness
    - kd_txt, ks_txt, kt_txt, rs_txt: texture indices (-1 for
    none)
    - use_phong: whether to use phong

### Function set_material_gltf_metallic_roughness()

~~~ .cpp
void set_material_gltf_metallic_roughness(scene* scn, int mid,
    const ym::vec3f& kb, float km, float rs, float op, int kd_txt, int km_txt);
~~~

Sets a gltf metallic roughness material reflectance.

- Parameters:
    - scn: scene
    - mid: material id
    - kb: base color term
    - km: metallic term
    - rs: specular roughness
    - kd_txt, km_txt: texture indices (-1 for none)

### Function set_material_gltf_specular_glossiness()

~~~ .cpp
void set_material_gltf_specular_glossiness(scene* scn, int mid,
    const ym::vec3f& kd, const ym::vec3f& ks, float rs, float op, int kd_txt,
    int ks_txt);
~~~

Sets a gltf metallic specular glossiness reflectance.

- Parameters:
    - scn: scene
    - mid: material id
    - ke: emission, term
    - kd: diffuse term
    - ks: specular term
    - rs: specular glossiness
    - kd_txt, ks_txt: texture indices (-1 for
    none)

### Function set_material_thin_glass()

~~~ .cpp
void set_material_thin_glass(scene* scn, int mid, const ym::vec3f& ks,
    const ym::vec3f& kt, int ks_txt, int kt_txt);
~~~

Sets a thin glass material

- Parameters:
    - scn: scene
    - mid: material id
    - ks: reflection
    - kt: transmission
    - ks_txt, kt_txt: texture indices (-1 for
    none)

### Function set_material_double_sided()

~~~ .cpp
void set_material_double_sided(scene* scn, int mid, bool double_sided);
~~~

Sets the material emission.

- Parameters:
    - scn: scene
    - mid: material id
    - double_sided: whether the material is double sided

### Function add_environment()

~~~ .cpp
int add_environment(
    scene* scn, const ym::frame3f& frame, const ym::vec3f& ke, int txt_id = -1);
~~~

Sets an environment in the scene.

- Parameters:
    - scn: scene
    - frame: local-to-world frame (x, y, z, o in column major order)
    - ke: emission
    - ke_txt: emission texture (-1 for none)
- Returns:
    - environment id

### Function add_triangle_shape()

~~~ .cpp
int add_triangle_shape(scene* scn, int ntriangles, const ym::vec3i* triangles,
    int nverts, const ym::vec3f* pos, const ym::vec3f* norm,
    const ym::vec2f* texcoord = nullptr, const ym::vec4f* color = nullptr,
    const ym::vec4f* tangsp = nullptr);
~~~

Sets a shape in the scene.

- Parameters:
    - scn: scene
    - frame: local-to-world frame (x, y, z, o in column major order)
    - mid: material id
    - npoints/nlines/ntriangles: number of elements
    - points/lines/tiangles: elem data
    - nverts: number of vertices
    - pos: vertex positions
    - norm/tang: vertex normals/tangents
    - texcoord: vertex texcoord
    - color: vertex color
    - tangsp: tangent space for normal and bump mapping
- Returns:
    - shape id

### Function add_point_shape()

~~~ .cpp
int add_point_shape(scene* scn, int npoints, const int* points, int nverts,
    const ym::vec3f* pos, const ym::vec3f* norm,
    const ym::vec2f* texcoord = nullptr, const ym::vec4f* color = nullptr,
    const float* radius = nullptr);
~~~

Sets a shape in the scene.

- Parameters:
    - scn: scene
    - sid: shape id
    - frame: local-to-world frame (x, y, z, o in column major order)
    - mid: material id
    - npoints/nlines/ntriangles: number of elements
    - points/lines/tiangles: elem data
    - nverts: number of vertices
    - pos: vertex positions
    - norm/tang: vertex normals/tangents
    - texcoord: vertex texcoord
    - color: vertex color
    - radius: vertex radius
- Returns:
    - shape id

### Function add_line_shape()

~~~ .cpp
int add_line_shape(scene* scn, int nlines, const ym::vec2i* lines, int nverts,
    const ym::vec3f* pos, const ym::vec3f* tang,
    const ym::vec2f* texcoord = nullptr, const ym::vec4f* color = nullptr,
    const float* radius = nullptr);
~~~

Sets a shape in the scene.

- Parameters:
    - scn: scene
    - frame: local-to-world frame (x, y, z, o in column major order)
    - mid: material id
    - npoints/nlines/ntriangles: number of elements
    - points/lines/tiangles: elem data
    - nverts: number of vertices
    - pos: vertex positions
    - norm/tang: vertex normals/tangents
    - texcoord: vertex texcoord
    - color: vertex color
    - radius: vertex radius
- Returns:
    - shape id

### Function add_instance()

~~~ .cpp
int add_instance(scene* scn, const ym::frame3f& frame, int sid, int mid);
~~~

Adds an instance in the scene.

- Parameters:
    - scn: scene
    - iid: instance id
    - frame: local-to-world frame (x, y, z, o in column major order)
    - sid: shape id
    - mid: material id
- Returns:
    - instance id

### Function specular_exponent_to_roughness()

~~~ .cpp
float specular_exponent_to_roughness(float n);
~~~

Convert a Phong exponent to GGX/Phong roughness

### Function specular_fresnel_from_ks()

~~~ .cpp
void specular_fresnel_from_ks(
    const ym::vec3f& ks, ym::vec3f& es, ym::vec3f& esk);
~~~

Estimates the fresnel coefficient es from ks at normal incidence

### Struct intersect_point

~~~ .cpp
struct intersect_point {
    float dist = 0;
    int iid = -1;
    int sid = -1;
    int eid = -1;
    ym::vec3f euv = {0, 0, 0};
    operator bool() const; 
}
~~~

Ray-scene Intersection.

- Members:
    - dist:      ray distance
    - iid:      instance index
    - sid:      shape index
    - eid:      element distance
    - euv:      element baricentric coordinates
    - operator bool():      check whether it was a hit


### Typedef intersect_first_cb

~~~ .cpp
using intersect_first_cb = std::function<intersect_point(const ym::ray3f& ray)>;
~~~

Ray-scene closest intersection callback.

- Parameters:
    - ray: ray
- Return:
    - intersection point

### Typedef intersect_any_cb

~~~ .cpp
using intersect_any_cb = std::function<bool(const ym::ray3f& ray)>;
~~~

Ray-scene intersection callback

- Parameters:
    - ray: ray
- Return:
    - whether we intersect or not

### Function set_intersection_callbacks()

~~~ .cpp
void set_intersection_callbacks(scene* scn, void* ctx,
    intersect_first_cb intersect_first, intersect_any_cb intersect_any);
~~~

Sets the intersection callbacks

### Function init_intersection()

~~~ .cpp
void init_intersection(scene* scn);
~~~

Initialize acceleration structure.

- Parameters:
    - scn: trace scene

### Typedef logging_cb

~~~ .cpp
using logging_cb = std::function<void(const char*)>;
~~~

Logging callback

### Function set_logging_callbacks()

~~~ .cpp
void set_logging_callbacks(
    scene* scn, logging_cb log_info, logging_cb log_error);
~~~

Logger

### Function init_lights()

~~~ .cpp
void init_lights(scene* scn);
~~~

Initialize lighting.

- Parameters:
    - scn: trace scene

### Enum shader_type

~~~ .cpp
enum struct shader_type {
    pathtrace = 0,
    eyelight,
    direct,
    direct_ao,
}
~~~

Type of rendering algorithm (shader)

- Values:
    - pathtrace:      pathtrace
    - eyelight:      eye hight for quick previews
    - direct:      direct illumination
    - direct_ao:      direct illumination with ambient occlusion


### Enum rng_type

~~~ .cpp
enum struct rng_type {
    uniform = 0,
    stratified,
    cmjs,
}
~~~

Random number generator type

- Values:
    - uniform:      uniform random numbers
    - stratified:      stratified random numbers
    - cmjs:      correlated multi-jittered sampling


### Enum filter_type

~~~ .cpp
enum struct filter_type {
    box = 1,
    triangle = 2,
    cubic = 3,
    catmull_rom = 4,
    mitchell = 5

}
~~~

Filter type

- Values:
    - box:      box filter
    - triangle:      hat filter
    - cubic:      cubic spline
    - catmull_rom:      Catmull-Rom spline
    - mitchell:      Mitchell-Netrevalli


### Struct trace_params

~~~ .cpp
struct trace_params {
    int camera_id = 0;
    int width = 0;
    int height = 0;
    int nsamples = 256;
    shader_type stype = shader_type::pathtrace;
    rng_type rtype = rng_type::stratified;
    filter_type ftype = filter_type::box;
    bool aux_buffers = false;
    ym::vec3f amb = {0, 0, 0};
    bool envmap_invisible = false;
    int min_depth = 3;
    int max_depth = 8;
    float pixel_clamp = 10;
    float ray_eps = 1e-4f;
    bool parallel = true;
}
~~~

Rendering params

- Members:
    - camera_id:      camera id
    - width:      width
    - height:      height
    - nsamples:      number of samples
    - stype:      sampler type
    - rtype:      random number generation type
    - ftype:      filter type
    - aux_buffers:      compute auxiliary buffers
    - amb:      ambient lighting
    - envmap_invisible:      view environment map
    - min_depth:      minimum ray depth
    - max_depth:      maximum ray depth
    - pixel_clamp:      final pixel clamping
    - ray_eps:      ray intersection epsilon
    - parallel:      parallel execution


### Function trace_block()

~~~ .cpp
void trace_block(const scene* scn, ym::vec4f* img, int block_x, int block_y,
    int block_width, int block_height, int samples_min, int samples_max,
    const trace_params& params);
~~~

Renders a block of sample

Notes: It is safe to call the function in parallel on different blocks.
But two threads should not access the same pixels at the same time.
Also blocks with different samples should be called sequentially if
accumulate is true.

- Parameters:
    - scn: trace scene
    - cid: camera id
    - img: pixel data in RGBA format (width/height in params)
    - nsamples: number of samples
    - block_x, block_y: block corner
    - block_width, block_height: block width and height
    - samples_min, samples_max: sample block to render
      [sample_min, sample_max]; max values are excluded

### Struct trace_state

~~~ .cpp
struct trace_state;
~~~

State for progressive rendering and denoising

### Function make_state()

~~~ .cpp
trace_state* make_state();
~~~

Creates state

### Function init_state()

~~~ .cpp
void init_state(
    trace_state* state, const scene* scn, const trace_params& params);
~~~

Initialize state

### Function free_state()

~~~ .cpp
void free_state(trace_state*& state);
~~~

Clear state

### Function get_traced_image()

~~~ .cpp
ym::image4f& get_traced_image(trace_state* state);
~~~

Grabs a reference to the image from the state

### Function get_aux_buffers()

~~~ .cpp
void get_aux_buffers(const trace_state* state, ym::image4f& norm,
    ym::image4f& albedo, ym::image4f& depth);
~~~

Grabs the image from the state

### Function get_cur_sample()

~~~ .cpp
int get_cur_sample(const trace_state* state);
~~~

Gets the current sample number

### Function trace_next_samples()

~~~ .cpp
bool trace_next_samples(trace_state* state, int nsamples);
~~~

Trace the next nsamples samples.

### Function trace_image()

~~~ .cpp
inline ym::image4f trace_image(const scene* scn, const trace_params& params);
~~~

Trace the whole image

### Function trace_async_start()

~~~ .cpp
void trace_async_start(trace_state* state);
~~~

Starts an anyncrhounous renderer with a maximum of 256 samples.

### Function trace_async_stop()

~~~ .cpp
void trace_async_stop(trace_state* state);
~~~

Stop the asynchronous renderer.

