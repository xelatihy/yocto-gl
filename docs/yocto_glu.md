# Yocto/Glu

A set of utilities to draw on screen with OpenGL 3.3. Mostly
used to make quick viewers. Not sure it is helpful to others (OpenGL is
really a portability mess, but there is nothing else).

This library depends in yocto_math.h
The library depends on GLEW for OpenGL functions on Windows and Linux.
Soon porting to glad.


## Features

1. image viewing with `shade_image()`, with support for tone mapping.
2. texture utilies to quickly create/update textures
    - create textures with `make_texture()`
    - create textures with `update_texture()`
    - delete textures with `clear_texture()`
3. program utilities in modern namespace
    - buffer objects with `create_buffer()`, `update_buffer()`
    - program creation/cleaning with `make_program()`, `clear_program()`
    - uniforms with `set_uniform()`
    - vertex attrib with `set_vertattr()` and `set_vertattr_buff()`
    - `draw_elems()`
4. a standard shader for GGX fragment shading and multiple lights in
  the `stdshader` namespace
    - initialize the shaders with `make_program()`
    - start/end each frame with `begin_frame()`, `end_frame()`
    - define lights with `set_lights()`
    - start/end each shape with `begin_shape()`, `end_shape()`
    - define material Parameters with `set_material()`
    - define vertices with `set_vert()`
    - draw elements with `draw_elems()`

The interface for each function is described in details in the interface
section of this file.

## History

- v 0.15: reference interface for images
- v 0.14: alpha cut out in stdshader
- v 0.13: simpler texture creation functions
- v 0.12: removing legacy functions
- v 0.11: use yocto_math in the interface and remove inline compilation
- v 0.10: user interface with dear ImGui
- v 0.9: user interface with GLFW and NUKLEAR
- v 0.8: switch to .h/.cpp pair
- v 0.7: cleaner srgb support for stdshader
- v 0.6: doxygen comments
- v 0.5: support for OpenGL 3.2
- v 0.4: support for legacy OpenGL
- v 0.3: [major API change] move to modern C++ interface
- v 0.2: removal of C interface
- v 0.1: C/C++ implementation
- v 0.0: initial release in C99

## Namespace yglu

OpenGL abstraction

### Typedef uint

~~~ .cpp
using uint = unsigned int;
~~~

Shortcut for GLuint.

### Typedef byte

~~~ .cpp
using byte = unsigned char;
~~~

Shortcut for unsigned chars.

### Enum etype : int

~~~ .cpp
enum struct etype : int {
    point = 1,
    line = 2,
    triangle = 3,
    quad = 4,
}
~~~

Shape types

- Values:
    - point:      points
    - line:      lines
    - triangle:      triangles
    - quad:      quads


### Enum ltype : int

~~~ .cpp
enum struct ltype : int {
    point = 0,
    directional = 1,
}
~~~

Light types

- Values:
    - point:      point lights
    - directional:      directional lights


### Enum texture_wrap

~~~ .cpp
enum struct texture_wrap {
    not_set = 0,
    repeat = 1,
    clamp = 2,
    mirror = 3,
}
~~~

Wrap values for texture

- Values:
    - not_set:      not set
    - repeat:      repeat
    - clamp:      clamp to edge
    - mirror:      mirror


### Enum texture_filter

~~~ .cpp
enum struct texture_filter {
    not_set = 0,
    linear = 1,
    nearest = 2,
    linear_mipmap_linear = 3,
    nearest_mipmap_nearest = 4,
    linear_mipmap_nearest = 5,
    nearest_mipmap_linear = 6,
}
~~~

Filter values for texture

- Values:
    - not_set:      not set
    - linear:      linear
    - nearest:      nearest
    - linear_mipmap_linear:      mip-mapping
    - nearest_mipmap_nearest:      mip-mapping
    - linear_mipmap_nearest:      mip-mapping
    - nearest_mipmap_linear:      mip-mapping


### Struct texture_info

~~~ .cpp
struct texture_info {
    uint txt_id = 0;
    int texcoord = 0;
    float scale = 1;
    texture_wrap wrap_s = texture_wrap::not_set;
    texture_wrap wrap_t = texture_wrap::not_set;
    texture_filter filter_mag = texture_filter::not_set;
    texture_filter filter_min = texture_filter::not_set;
    texture_info(); 
    texture_info(uint tid); 
}
~~~

Texture information for parameter setting.

- Members:
    - txt_id:      texture id
    - texcoord:      texture coordinate set
    - scale:      texture strength/scale (used by some models)
    - wrap_s:      wrap mode
    - wrap_t:      wrap mode
    - filter_mag:      filter mode
    - filter_min:      filter mode
    - texture_info():      default constructor
    - texture_info():      constructor from texture id only


### Enum buffer_type

~~~ .cpp
enum struct buffer_type {
    vertex = 0,
    element = 1,
    uniform_block = 2

}
~~~

Buffer type

- Values:
    - vertex:      vertex
    - element:      element
    - uniform_block:      uniform block


### Function check_error()

~~~ .cpp
bool check_error(bool print = true);
~~~

Checks for GL error and then prints

### Function clear_buffers()

~~~ .cpp
void clear_buffers(const ym::vec4f& background =;
~~~

Clear window

### Function enable_depth_test()

~~~ .cpp
void enable_depth_test(bool enabled);
~~~

Enable/disable depth test

### Function enable_culling()

~~~ .cpp
void enable_culling(bool enabled);
~~~

Enable/disable culling

### Function enable_wireframe()

~~~ .cpp
void enable_wireframe(bool enabled);
~~~

Enable/disable wireframe

### Function enable_edges()

~~~ .cpp
void enable_edges(bool enabled, float tolerance = 0.9999f);
~~~

Enable/disable edges. Attempts to avoid z-fighting but the method is not
robust.

### Function enable_blending()

~~~ .cpp
void enable_blending(bool enabled);
~~~

Enable/disable blending

### Function set_blend_over()

~~~ .cpp
void set_blend_over();
~~~

Set blending to over operator

### Function line_width()

~~~ .cpp
void line_width(float w);
~~~

Line width

### Function set_viewport()

~~~ .cpp
void set_viewport(const ym::vec4i& v);
~~~

Set viewport

### Function shade_image()

~~~ .cpp
void shade_image(
    uint tid, int win_w, int win_h, float ox, float oy, float zoom);
~~~

Draw an texture tid of size img_w, img_h on a window of size win_w, win_h
with top-left corner at ox, oy with a zoom zoom.

### Function shade_image()

~~~ .cpp
void shade_image(uint tid, int win_w, int win_h, float ox, float oy, float zoom,
    ym::tonemap_type tmtype, float exposure, float gamma);
~~~

As above but includes an exposure/gamma correction.

### Function make_texture()

~~~ .cpp
uint make_texture(int w, int h, int nc, const float* pixels, bool linear,
    bool mipmap, bool as_float);
~~~

Creates a texture with pixels values of size w, h with nc number of
components (1-4).
Internally use float if as_float and filtering if filter.
Returns the texture id.

### Function make_texture()

~~~ .cpp
uint make_texture(int w, int h, int nc, const unsigned char* pixels,
    bool linear, bool mipmap, bool as_srgb);
~~~

Creates a texture with pixels values of size w, h with nc number of
components (1-4).
Internally use srgb lookup if as_srgb and filtering if filter.
Returns the texture id.

### Function update_texture()

~~~ .cpp
void update_texture(
    uint tid, int w, int h, int nc, const float* pixels, bool mipmap);
~~~

Updates the texture tid with new image data.

### Function update_texture()

~~~ .cpp
void update_texture(
    uint tid, int w, int h, int nc, const unsigned char* pixels, bool mipmap);
~~~

Updates the texture tid with new image data.

### Function make_texture()

~~~ .cpp
template <int N>
inline uint make_texture(const ym::image<ym::vec<float, N>>& img, bool linear,
    bool mipmap, bool as_float);
~~~

Creates a texture from an image.
Internally use float if as_float and filtering if filter.
Returns the texture id.

### Function make_texture()

~~~ .cpp
template <int N>
inline uint make_texture(const ym::image<ym::vec<byte, N>>& img, bool linear,
    bool mipmap, bool as_srgb);
~~~

Creates a texture from an image.
Internally use srgb lookup if as_srgb and filtering if filter.
Returns the texture id.

### Function update_texture()

~~~ .cpp
template <int N>
inline void update_texture(
    uint tid, const ym::image<ym::vec<float, N>>& img, bool mipmap);
~~~

Updates the texture tid with new image data.

### Function update_texture()

~~~ .cpp
template <int N>
inline void update_texture(
    uint tid, const ym::image<ym::vec<byte, N>>& img, bool mipmap);
~~~

Updates the texture tid with new image data.

### Function clear_texture()

~~~ .cpp
void clear_texture(uint* tid);
~~~

Destroys the texture tid.

### Function make_buffer_raw()

~~~ .cpp
uint make_buffer_raw(
    int size, const void* values, buffer_type btype, bool dynamic);
~~~

Creates a buffer with num elements of size size stored in values, where
content is dyanamic if dynamic.
Returns the buffer id.

### Function make_buffer()

~~~ .cpp
template <typename T>
inline uint make_buffer(int num, const T* values, bool elements, bool dynamic);
~~~

Creates a buffer with num elements of size size stored in values, where
content is dyanamic if dynamic.
Returns the buffer id.

### Function make_buffer()

~~~ .cpp
template <typename T>
inline uint make_buffer(
    const std::vector<T>& values, bool elements, bool dynamic);
~~~

Creates a buffer with num elements of size size stored in values, where
content is dyanamic if dynamic.
Returns the buffer id.

### Function make_buffer()

~~~ .cpp
template <typename T>
inline uint make_buffer(const T& values, buffer_type btype, bool dynamic);
~~~

Creates a buffer with a block of data stored in values, where
content is dyanamic if dynamic.
Returns the buffer id.

### Function update_buffer_raw()

~~~ .cpp
void update_buffer_raw(
    uint bid, int size, const void* values, buffer_type btype, bool dynamic);
~~~

Updates the buffer bid with new data.

### Function update_buffer()

~~~ .cpp
template <typename T>
inline void update_buffer(
    uint bid, int num, const T* values, bool elements, bool dynamic);
~~~

Updates the buffer bid with new data.

### Function update_buffer()

~~~ .cpp
template <typename T>
inline void update_buffer(
    uint bid, const std::vector<T>& values, bool elements, bool dynamic);
~~~

Updates the buffer bid with new data.

### Function update_buffer()

~~~ .cpp
template <typename T>
inline void update_buffer(
    uint bid, const T& values, buffer_type btype, bool dynamic);
~~~

Updates the buffer bid with new data.

### Function clear_buffer()

~~~ .cpp
void clear_buffer(uint* bid);
~~~

Destroys the buffer bid.

### Function make_vertex_arrays()

~~~ .cpp
uint make_vertex_arrays();
~~~

Creates and OpenGL vertex array object.

### Function clear_vertex_arrays()

~~~ .cpp
void clear_vertex_arrays(uint* aid);
~~~

Destroys the program pid and optionally the sahders vid and fid.

### Function bind_vertex_arrays()

~~~ .cpp
void bind_vertex_arrays(uint aid);
~~~

Binds a vertex array object

### Function make_program()

~~~ .cpp
uint make_program(const std::string& vertex, const std::string& fragment,
    uint* vid, uint* fid, uint* vao);
~~~

Creates and OpenGL program from vertex and fragment code. Returns the
program id. Optionally return vertex and fragment shader ids. A VAO is
created.

### Function clear_program()

~~~ .cpp
void clear_program(uint* pid, uint* vid, uint* fid, uint* vao);
~~~

Destroys the program pid and optionally the sahders vid and fid.

### Function bind_program()

~~~ .cpp
void bind_program(uint pid);
~~~

Binds a program

### Function get_uniform_location()

~~~ .cpp
int get_uniform_location(uint pid, const std::string& name);
~~~

Get uniform location (simple GL wrapper that avoids GL includes)

### Function get_uniform_block_index()

~~~ .cpp
int get_uniform_block_index(uint pid, const std::string& name);
~~~

Get uniform location (simple GL wrapper that avoids GL includes)

### Function get_attrib_location()

~~~ .cpp
int get_attrib_location(uint pid, const std::string& name);
~~~

Get attrib location (simple GL wrapper that avoids GL includes)

### Function get_uniforms_names()

~~~ .cpp
std::vector<std::pair<std::string, int>> get_uniforms_names(uint pid);
~~~

Get the names of all uniforms

### Function get_uniforms_namemap()

~~~ .cpp
inline std::unordered_map<std::string, int> get_uniforms_namemap(uint pid);
~~~

Get the names of all uniforms

### Function get_uniform_buffers_names()

~~~ .cpp
std::vector<std::pair<std::string, int>> get_uniform_buffers_names(uint pid);
~~~

Get the names of all uniform blocks

### Function get_uniform_buffers_namemap()

~~~ .cpp
inline std::unordered_map<std::string, int> get_uniform_buffers_namemap(
    uint pid);
~~~

Get the names of all uniform blocks

### Function get_attributes_names()

~~~ .cpp
std::vector<std::pair<std::string, int>> get_attributes_names(uint pid);
~~~

Get the names of all attributes

### Function get_attributes_namemap()

~~~ .cpp
inline std::unordered_map<std::string, int> get_attributes_namemap(uint pid);
~~~

Get the names of all attributes

### Function set_uniform_raw()

~~~ .cpp
bool set_uniform_raw(uint prog, int var, const int* val, int nc, int count);
~~~

Set uniform integer values val for program prog and variable loc.
The values have nc number of components (1-4) and count elements
(for arrays).

### Function set_uniform_raw()

~~~ .cpp
bool set_uniform_raw(uint prog, int var, const float* val, int nc, int count);
~~~

Set uniform float values val for program prog and variable var.
The values have nc number of components (1-4) and count elements
(for arrays).

### Function set_uniform()

~~~ .cpp
inline bool set_uniform(uint prog, int var, bool val);
~~~

Set uniform float values val for program prog and variable var.

### Function set_uniform()

~~~ .cpp
inline bool set_uniform(uint prog, int var, int val);
~~~

Set uniform float values val for program prog and variable var.

### Function set_uniform()

~~~ .cpp
inline bool set_uniform(uint prog, int var, float val);
~~~

Set uniform float values val for program prog and variable var.

### Function set_uniform()

~~~ .cpp
template <int N>
inline bool set_uniform(uint prog, int var, const ym::vec<float, N>& val);
~~~

Set uniform float values val for program prog and variable var.

### Function set_uniform()

~~~ .cpp
template <int N>
inline bool set_uniform(uint prog, int var, const ym::vec<int, N>& val);
~~~

Set uniform float values val for program prog and variable var.

### Function set_uniform()

~~~ .cpp
inline bool set_uniform(uint prog, int var, const ym::mat4f& val);
~~~

Set uniform float values val for program prog and variable var.

### Function set_uniform()

~~~ .cpp
inline bool set_uniform(uint prog, int var, const ym::frame3f& val);
~~~

Set uniform float values val for program prog and variable var.

### Function set_uniform()

~~~ .cpp
inline bool set_uniform(uint prog, int var, const int* val, int num);
~~~

Set uniform float values val for program prog and variable var.

### Function set_uniform()

~~~ .cpp
inline bool set_uniform(uint prog, int var, const float* val, int num);
~~~

Set uniform float values val for program prog and variable var.

### Function set_uniform()

~~~ .cpp
template <int N>
inline bool set_uniform(
    uint prog, int var, const ym::vec<float, N>* val, int num);
~~~

Set uniform float values val for program prog and variable var.

### Function set_uniform()

~~~ .cpp
template <int N>
inline bool set_uniform(
    uint prog, int var, const ym::vec<int, N>* val, int num);
~~~

Set uniform float values val for program prog and variable var.

### Function set_uniform()

~~~ .cpp
inline bool set_uniform(uint prog, int var, const ym::mat4f* val, int num);
~~~

Set uniform float values val for program prog and variable var.

### Function set_uniform()

~~~ .cpp
template <typename T>
inline bool set_uniform(uint prog, const std::string& var, const T& val);
~~~

Set uniform float values val for program prog and variable var.

### Function set_uniform()

~~~ .cpp
template <typename T>
inline bool set_uniform(
    uint prog, const std::string& var, const T* val, int num);
~~~

Set uniform float values val for program prog and variable var.

### Function set_uniform_block()

~~~ .cpp
bool set_uniform_block(uint prog, int var, int bid, uint bunit);
~~~

Set uniform block bid for program prog and variable var.

### Function set_uniform_texture()

~~~ .cpp
bool set_uniform_texture(
    uint prog, int var, const texture_info& tinfo, uint tunit);
~~~

Set uniform texture id tid and unit tunit for program prog and variable var.

### Function set_uniform_texture()

~~~ .cpp
inline bool set_uniform_texture(
    uint prog, int var, int varon, const texture_info& tinfo, uint tunit);
~~~

Set uniform texture id tid and unit tunit for program prog and variable var.
Optionally sets the int variable varon to 0/1 whether the texture is enable
on not.

### Function set_uniform_texture()

~~~ .cpp
inline bool set_uniform_texture(
    uint prog, const std::string& var, const texture_info& tinfo, uint tunit);
~~~

Set uniform texture id tid and unit tunit for program prog and variable var.

### Function set_uniform_texture()

~~~ .cpp
inline bool set_uniform_texture(uint prog, const std::string& var,
    const std::string& varon, const texture_info& tinfo, uint tunit);
~~~

Set uniform texture id tid and unit tunit for program prog and variable var.
Optionally sets the int variable varon to 0/1 whether the texture is enable
on not.

### Function set_vertattr_val_raw()

~~~ .cpp
bool set_vertattr_val_raw(uint prog, int var, const float* value, int nc);
~~~

Sets a constant value for a vertex attribute for program prog and
variable var. The attribute has nc components.

### Function set_vertattr_val_raw()

~~~ .cpp
bool set_vertattr_val_raw(uint prog, int var, const int* value, int nc);
~~~

Sets a constant value for a vertex attribute for program prog and
variable var. The attribute has nc components.

### Function set_vertattr_buffer()

~~~ .cpp
bool set_vertattr_buffer(uint prog, int var, uint bid, int nc);
~~~

Sets a vartex attribute for program prog and variable var to the buffer bid.
The attribute has nc components and per-vertex values values.

### Function set_vertattri_buffer()

~~~ .cpp
bool set_vertattri_buffer(uint prog, int var, uint bid, int nc);
~~~

Sets a vartex attribute for program prog and variable var to the buffer bid.
The attribute has nc components and per-vertex values values.

### Function set_vertattr_raw()

~~~ .cpp
bool set_vertattr_raw(uint prog, int var, uint bid, int nc, const float* def);
~~~

Sets a vartex attribute for program prog and variable var. The attribute
has nc components and either buffer bid or a single value def
(if bid is zero). Convenience wrapper to above functions.

### Function set_vertattr_raw()

~~~ .cpp
bool set_vertattr_raw(uint prog, int var, uint bid, int nc, const int* def);
~~~

Sets a vartex attribute for program prog and variable var. The attribute
has nc components and either buffer bid or a single value def
(if bid is zero). Convenience wrapper to above functions.

### Function set_vertattr()

~~~ .cpp
template <int N>
inline bool set_vertattr(
    uint prog, int var, uint bid, const ym::vec<float, N>& def);
~~~

Sets a vartex attribute for program prog and variable var. The attribute
is either a buffer bid or a single value def
(if bid is zero). Convenience wrapper to above functions.

### Function set_vertattr()

~~~ .cpp
template <int N>
inline bool set_vertattr(
    uint prog, int var, uint bid, const ym::vec<int, N>& def);
~~~

Sets a vartex attribute for program prog and variable var. The attribute
is either a buffer bid or a single value def
(if bid is zero). Convenience wrapper to above functions.

### Function set_vertattr()

~~~ .cpp
template <int N>
inline bool set_vertattr(
    uint prog, const std::string& var, uint bid, const ym::vec<float, N>& def);
~~~

Sets a vartex attribute for program prog and variable var. The attribute
is either a buffer bid or a single value def
(if bid is zero). Convenience wrapper to above functions.

### Function set_vertattr()

~~~ .cpp
template <int N>
inline bool set_vertattr(
    uint prog, const std::string& var, uint bid, const ym::vec<int, N>& def);
~~~

Sets a vartex attribute for program prog and variable var. The attribute
is either a buffer bid or a single value def
(if bid is zero). Convenience wrapper to above functions.

### Function draw_elems()

~~~ .cpp
bool draw_elems(int nelems, uint bid, etype etype);
~~~

Draws nelems elements elem of type etype.

## Namespace stdshader

Shade with a physically-based standard shader based on Phong/GGX.

### Function make_program()

~~~ .cpp
uint make_program(uint* vao);
~~~

Initialize a standard shader.

### Function begin_frame()

~~~ .cpp
void begin_frame(uint prog, uint vao, bool shade_eyelight, float img_exposure,
    ym::tonemap_type img_tonemap, float img_gamma,
    const ym::mat4f& camera_xform, const ym::mat4f& camera_xform_inv,
    const ym::mat4f& camera_proj);
~~~

Starts a frame by setting exposure/gamma values, camera transforms and
projection. Sets also whether to use full shading or a quick eyelight
preview.

### Function end_frame()

~~~ .cpp
void end_frame();
~~~

Ends a frame.

### Function set_lights()

~~~ .cpp
void set_lights(uint prog, const ym::vec3f& amb, int num, ym::vec3f* pos,
    ym::vec3f* ke, ltype* ltype);
~~~

Set num lights with position pos, color ke, type ltype. Also set the
ambient illumination amb.

### Function begin_shape()

~~~ .cpp
void begin_shape(uint prog, const ym::mat4f& xform);
~~~

Begins drawing a shape with transform xform.

### Function end_shape()

~~~ .cpp
void end_shape();
~~~

End shade drawing.

### Function set_highlight()

~~~ .cpp
void set_highlight(uint prog, const ym::vec4f& highlight);
~~~

Set the object as highlight color.

### Function set_material_emission_only()

~~~ .cpp
void set_material_emission_only(uint prog, const ym::vec3f& ke, float op,
    const texture_info& ke_txt, bool double_sided, bool alpha_cutout);
~~~

Set material values for emission only (constant color).
Indicates textures ids with the correspoinding XXX_txt variables.
Works for points/lines/triangles. Element type set by draw_XXX calls.

### Function set_material_generic()

~~~ .cpp
void set_material_generic(uint prog, const ym::vec3f& ke, const ym::vec3f& kd,
    const ym::vec3f& ks, float rs, float op, const texture_info& ke_txt,
    const texture_info& kd_txt, const texture_info& ks_txt,
    const texture_info& rs_txt, const texture_info& norm_txt,
    const texture_info& occ_txt, bool use_phong, bool double_sided,
    bool alpha_cutout);
~~~

Set material values with emission ke, diffuse kd, specular ks and
specular roughness rs, opacity op. Indicates textures ids with the
correspoinding XXX_txt variables. Sets also normal and occlusion
maps. Works for points/lines/triangles (diffuse for points,
Kajiya-Kay for lines, GGX/Phong for triangles). Element type set by
draw_XXX calls.

### Function set_material_gltf_metallic_roughness()

~~~ .cpp
void set_material_gltf_metallic_roughness(uint prog, const ym::vec3f& ke,
    const ym::vec3f& kb, float km, float rs, float op,
    const texture_info& ke_txt, const texture_info& kb_txt,
    const texture_info& km_txt, const texture_info& norm_txt,
    const texture_info& occ_txt, bool use_phong, bool double_sided,
    bool alpha_cutout);
~~~

Set material values for glTF specular-roughness PBR shader,
with emission ke, base color kb, opacity op, metallicity km and
specular roughness rs. Uses basecolor-opacity texture kb_txt and
metallic-roughness texture km_txt. Sets also normal and occlusion
maps. Works for points/lines/triangles (diffuse for points, Kajiya-Kay
for lines, GGX/Phong for triangles). Element type set by draw_XXX calls.

### Function set_material_gltf_specular_glossiness()

~~~ .cpp
void set_material_gltf_specular_glossiness(uint prog, const ym::vec3f& ke,
    const ym::vec3f& kd, const ym::vec3f& ks, float rs, float op,
    const texture_info& ke_txt, const texture_info& kd_txt,
    const texture_info& ks_txt, const texture_info& norm_txt,
    const texture_info& occ_txt, bool use_phong, bool double_sided,
    bool alpha_cutout);
~~~

Set material values for glTF specular-roughness PBR shader,
with emission ke, diffuse color kd, opacity op, specular ks and
specular glossiness rs. Uses diffuse-opacity texture kd_txt and
specular-glpossiness texture ks_txt. Sets also normal and occlusion
maps. Works for points/lines/triangles (diffuse for points, Kajiya-Kay
for lines, GGX/Phong for triangles). Element type set by draw_XXX calls.

### Function specular_exponent_to_roughness()

~~~ .cpp
float specular_exponent_to_roughness(float n);
~~~

Convertes a phong exponent to roughness.

### Function set_vert()

~~~ .cpp
void set_vert(uint prog, const ym::vec3f* pos, const ym::vec3f* norm,
    const ym::vec2f* texcoord, const ym::vec4f* color, const ym::vec4f* tangsp);
~~~

Set vertex data with position pos, normals norm, texture coordinates
texcoord and per-vertex color color and tangent space tangsp.

### Function set_vert()

~~~ .cpp
void set_vert(
    uint prog, uint pos, uint norm, uint texcoord, uint color, uint tangsp);
~~~

Set vertex data with buffers for position pos, normals norm, texture
coordinates texcoord, per-vertex color color and tangent space tangsp.

### Function set_vert_skinning()

~~~ .cpp
void set_vert_skinning(
    uint prog, uint weights, uint joints, int nxforms, const ym::mat4f* xforms);
~~~

Set vertex data with buffers for skinning.

### Function set_vert_gltf_skinning()

~~~ .cpp
void set_vert_gltf_skinning(
    uint prog, uint weights, uint joints, int nxforms, const ym::mat4f* xforms);
~~~

Set vertex data with buffers for skinning.

### Function set_vert_skinning_off()

~~~ .cpp
void set_vert_skinning_off(uint prog);
~~~

Disables vertex skinning.

### Function draw_elems()

~~~ .cpp
void draw_elems(uint prog, int num, uint bid, etype etype);
~~~

Draw num elements elem of type etype.

### Function draw_points()

~~~ .cpp
void draw_points(uint prog, int num, uint bid);
~~~

Draw num elements elem of type etype.

### Function draw_lines()

~~~ .cpp
void draw_lines(uint prog, int num, uint bid);
~~~

Draw num elements elem of type etype.

### Function draw_triangles()

~~~ .cpp
void draw_triangles(uint prog, int num, uint bid);
~~~

Draw num elements elem of type etype.

