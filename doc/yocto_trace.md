# Yocto/Trace

Path tracer implementation for support for textured mesh
lights, GGX/Phong materials, environment mapping. The interface supports
progressive parallel execution with any splitting strategy, by
generating per-sample random number sequences or fully deterministic
hash-based sampling.

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
sequence per path. This means you can rul the path tracer in any order
serially or in parallel.

For now, we support a straightforward path tracer with explicit direct
illumination using MIS.

This library depends in yocto_math.h. Optionally depend on yocto_bvh.h/.cpp
for internal acceleration. Disable this by setting YTRACE_NO_BVH.


## Usage

1. define your scene setting up only the public data
    - init the scene with make_scene()
    - define cameras with set_camera()
    - defile shapes with set set_triagle_shape(), set_line_shape(),
set_point_shape()
    - define materials with set_material()
    - define textures with set_texture()
    - define environments with set_environment()
    - set intersection routines with set_intersection_callbacks()
  or use the builtin intersection with init_intersection()
2. prepare for rendering with init_lights()
3. define rendering params in render_params
4. render blocks of samples with trace_block() or the whole image with
    trace_image()


## History

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
void free_scene(scene* scn);
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
int add_texture(scene* scn, int width, int height, int ncomp, const float* hdr);
~~~

Adds a texture in the scene.

- Parameters:
    - scn: scene
    - tid: texture id
    - width: width
    - height: height
    - ncomp: number of components (1-4)
    - hdr: hdr pixels
    - ldr: ldr pixels (sRGB)
- Returns:
    - texture id

### Function add_texture()

~~~ .cpp
int add_texture(scene* scn, int width, int height, int ncomp, const byte* ldr);
~~~

Sets a texture in the scene.

- Parameters:
    - scn: scene
    - tid: texture id
    - width: width
    - height: height
    - ncomp: number of components (1-4)
    - hdr: hdr pixels
    - ldr: ldr pixels (sRGB)
- Returns:
    - texture id

### Function add_material()

~~~ .cpp
int add_material(scene* scn, const ym::vec3f& ke, const ym::vec3f& kd,
    const ym::vec3f& ks, const ym::vec3f& kt, float rs = 0.1, int ke_txt = -1,
    int kd_txt = -1, int ks_txt = -1, int kt_txt = -1, int rs_txt = -1,
    int norm_txt = -1, bool use_phong = false);
~~~

Sets a material in the scene. [DEPRECATED]

- Parameters:
    - scn: scene
    - ke: emission, term
    - kd: diffuse term
    - ks: specular term
    - rs: specular roughness
    - ke_txt, kd_txt, ks_txt, rs_txt, norm_txt: texture indices (-1 for
    none) - use_phong: whether to use phong
- Returns:
    - material id

### Function add_material_generic()

~~~ .cpp
int add_material_generic(scene* scn, const ym::vec3f& ke, const ym::vec3f& kd,
    const ym::vec3f& ks, const ym::vec3f& kt, float rs, float op, int ke_txt,
    int kd_txt, int ks_txt, int kt_txt, int rs_txt, int op_txt, int norm_txt,
    int occ_txt, bool use_phong = false);
~~~

Sets a material in the scene with the most customization possible.

- Parameters:
    - scn: scene
    - ke: emission, term
    - kd: diffuse term
    - ks: specular term
    - rs: specular roughness
    - ke_txt, kd_txt, ks_txt, rs_txt, norm_txt: texture indices (-1 for
    none) - use_phong: whether to use phong
- Returns:
    - material id

### Function add_material_gltf_metallic_roughness()

~~~ .cpp
int add_material_gltf_metallic_roughness(scene* scn, const ym::vec3f& ke,
    const ym::vec3f& kd, float ks, float rs, float op, int ke_txt, int kd_txt,
    int ks_txt, int norm_txt, int occ_txt, bool use_phong = false);
~~~

Sets a gltf metallic roughness material.

- Parameters:
    - scn: scene
    - ke: emission, term
    - kd: diffuse term
    - ks: specular term
    - rs: specular roughness
    - ke_txt, kd_txt, ks_txt, rs_txt, norm_txt: texture indices (-1 for
    none)
    - use_phong: whether to use phong
- Returns:
    - material id

### Function add_material_gltf_specular_glossiness()

~~~ .cpp
int add_material_gltf_specular_glossiness(scene* scn, const ym::vec3f& ke,
    const ym::vec3f& kd, const ym::vec3f& ks, float rs, float op, int ke_txt,
    int kd_txt, int ks_txt, int norm_txt, int occ_txt, bool use_phong = false);
~~~

Sets a gltf metallic specular glossiness.

- Parameters:
    - scn: scene
    - ke: emission, term
    - kd: diffuse term
    - ks: specular term
    - rs: specular roughness
    - ke_txt, kd_txt, ks_txt, rs_txt, norm_txt: texture indices (-1 for
    none)
   - use_phong: whether to use phong
- Returns:
    - material id

### Function add_material_emission_only()

~~~ .cpp
int add_material_emission_only(
    scene* scn, const ym::vec3f& ke, int ke_txt, int norm_txt, int occ_txt);
~~~

Sets a gltf emission only material.

- Parameters:
    - scn: scene
    - mid: material id
    - ke: emission, term
    - kd: diffuse term
    - ks: specular term
    - rs: specular roughness
    - ke_txt, kd_txt, ks_txt, rs_txt, norm_txt: texture indices (-1 for
    none)
    - use_phong: whether to use phong
- Returns:
    - material id

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

### Function set_vert_material()

~~~ .cpp
void set_vert_material(scene* scn, int sid, const ym::vec3f* ke,
    const ym::vec3f* kd, const ym::vec3f* ks, const float* rs);
~~~

Sets per-vertex material properties.

- Parameters:
    - scn: scene
    - sid: shape id
    - ke: per-vertex emission
    - kd: per-vertex diffuse
    - ks: per-vertex specular
    - rs: per-vertex roughness

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
using intersect_first_cb = intersect_point (*)(void* ctx, const ym::ray3f& ray);
~~~

Ray-scene closest intersection callback.

- Parameters:
    - ctx: context
    - ray: ray
- Return:
    - intersection point

### Typedef intersect_any_cb

~~~ .cpp
using intersect_any_cb = bool (*)(void* ctx, const ym::ray3f& o);
~~~

Ray-scene intersection callback

- Parameters:
    - ctx: context
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

### Typedef logging_msg_cb

~~~ .cpp
using logging_msg_cb = void (*)(
    int level, const char* name, const char* msg, va_list args);
~~~

Logger callback

### Function set_logging_callbacks()

~~~ .cpp
void set_logging_callbacks(scene* scn, void* ctx, logging_msg_cb log);
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


### Struct render_params

~~~ .cpp
struct render_params {
    int camera_id = 0;
    int nsamples = 256;
    bool progressive = true;
    shader_type stype = shader_type::pathtrace;
    rng_type rtype = rng_type::stratified;
    ym::vec3f amb = {0, 0, 0};
    bool envmap_invisible = false;
    int min_depth = 3;
    int max_depth = 8;
    float pixel_clamp = 10;
    float ray_eps = 1e-4f;
}
~~~

Rendering params

- Members:
    - camera_id:      camera id
    - nsamples:      number of samples
    - progressive:      progressive rendering
    - stype:      sampler type
    - rtype:      random number generation type
    - amb:      ambient lighting
    - envmap_invisible:      view environment map
    - min_depth:      minimum ray depth
    - max_depth:      maximum ray depth
    - pixel_clamp:      final pixel clamping
    - ray_eps:      ray intersection epsilon


### Function trace_block()

~~~ .cpp
void trace_block(const scene* scn, int width, int height, ym::vec4f* img,
    int block_x, int block_y, int block_width, int block_height,
    int samples_min, int samples_max, const render_params& params);
~~~

Renders a block of sample

Notes: It is safe to call the function in parallel one different blocks.
But two threads should not access the same pixels at the same time.
Also blocks with different samples should be called sequentially if
accumulate is true.

- Parameters:
    - scn: trace scene
    - cid: camera id
    - width: image width
    - height: image height
    - img: pixel data in RGBA format
    - nsamples: number of samples
    - block_x, block_y: block corner
    - block_width, block_height: block width and height
    - samples_min, samples_max: sample block to render
      [sample_min, sample_max]; max values are excluded

### Function trace_image()

~~~ .cpp
void trace_image(const scene* scn, const int width, int height, ym::vec4f* img,
    const render_params& params);
~~~

Convenience function to call trace_block with all sample at once.

