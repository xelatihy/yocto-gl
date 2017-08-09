///
/// # Yocto/Glu
///
/// A set of utilities to draw on screen with OpenGL 3.3. Mostly
/// used to make quick viewers. Not sure it is helpful to others (OpenGL is
/// really a portability mess, but there is nothing else).
///
/// This library depends in yocto_math.h
/// The library depends on GLEW for OpenGL functions on Windows and Linux.
/// Soon porting to glad.
///
///
/// ## Features
///
/// 1. image viewing with `shade_image()`, with support for tone mapping.
/// 2. texture utilies to quickly create/update textures
///     - create textures with `make_texture()`
///     - create textures with `update_texture()`
///     - delete textures with `clear_texture()`
/// 3. program utilities in modern namespace
///     - buffer objects with `create_buffer()`, `update_buffer()`
///     - program creation/cleaning with `make_program()`, `clear_program()`
///     - uniforms with `set_uniform()`
///     - vertex attrib with `set_vertattr()` and `set_vertattr_buff()`
///     - `draw_elems()`
/// 4. a standard shader for GGX fragment shading and multiple lights in
///   the `stdshader` namespace
///     - initialize the shaders with `make_program()`
///     - start/end each frame with `begin_frame()`, `end_frame()`
///     - define lights with `set_lights()`
///     - start/end each shape with `begin_shape()`, `end_shape()`
///     - define material Parameters with `set_material()`
///     - define vertices with `set_vert()`
///     - draw elements with `draw_elems()`
///
/// The interface for each function is described in details in the interface
/// section of this file.
///
/// ## History
///
/// - v 0.15: reference interface for images
/// - v 0.14: alpha cut out in stdshader
/// - v 0.13: simpler texture creation functions
/// - v 0.12: removing legacy functions
/// - v 0.11: use yocto_math in the interface and remove inline compilation
/// - v 0.10: user interface with dear ImGui
/// - v 0.9: user interface with GLFW and NUKLEAR
/// - v 0.8: switch to .h/.cpp pair
/// - v 0.7: cleaner srgb support for stdshader
/// - v 0.6: doxygen comments
/// - v 0.5: support for OpenGL 3.2
/// - v 0.4: support for legacy OpenGL
/// - v 0.3: [major API change] move to modern C++ interface
/// - v 0.2: removal of C interface
/// - v 0.1: C/C++ implementation
/// - v 0.0: initial release in C99
///
namespace yglu {}

//
// LICENSE:
//
// Copyright (c) 2016 -- 2017 Fabio Pellacini
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//

#ifndef _YGLU_H_
#define _YGLU_H_

#include <array>
#include <string>
#include <unordered_map>
#include <vector>

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// INTERFACE
// -----------------------------------------------------------------------------

///
/// OpenGL abstraction
///
namespace yglu {

///
/// Shortcut for GLuint.
///
using uint = unsigned int;

///
/// Shortcut for unsigned chars.
///
using byte = unsigned char;

///
/// Shape types
///
enum struct etype : int {
    /// points
    point = 1,
    /// lines
    line = 2,
    /// triangles
    triangle = 3,
    /// quads
    quad = 4,
};

///
/// Light types
///
enum struct ltype : int {
    /// point lights
    point = 0,
    /// directional lights
    directional = 1,
};

///
/// Wrap values for texture
///
enum struct texture_wrap {
    /// not set
    not_set = 0,
    /// repeat
    repeat = 1,
    /// clamp to edge
    clamp = 2,
    /// mirror
    mirror = 3,
};

///
/// Filter values for texture
///
enum struct texture_filter {
    /// not set
    not_set = 0,
    /// linear
    linear = 1,
    /// nearest
    nearest = 2,
    /// mip-mapping
    linear_mipmap_linear = 3,
    /// mip-mapping
    nearest_mipmap_nearest = 4,
    /// mip-mapping
    linear_mipmap_nearest = 5,
    /// mip-mapping
    nearest_mipmap_linear = 6,
};

///
/// Texture information for parameter setting.
///
struct texture_info {
    /// texture id
    uint txt_id = 0;
    /// texture coordinate set
    int texcoord = 0;
    /// texture strength/scale (used by some models)
    float scale = 1;
    /// wrap mode
    texture_wrap wrap_s = texture_wrap::not_set;
    /// wrap mode
    texture_wrap wrap_t = texture_wrap::not_set;
    /// filter mode
    texture_filter filter_mag = texture_filter::not_set;
    /// filter mode
    texture_filter filter_min = texture_filter::not_set;

    /// default constructor
    texture_info() {}
    /// constructor from texture id only
    texture_info(uint tid) : txt_id(tid) {}
};

///
/// Buffer type
///
enum struct buffer_type {
    /// vertex
    vertex = 0,
    /// element
    element = 1,
    /// uniform block
    uniform_block = 2
};

///
/// Checks for GL error and then prints
///
bool check_error(bool print = true);

///
/// Clear window
///
void clear_buffers(const ym::vec4f& background = {0, 0, 0, 0});

///
/// Enable/disable depth test
///
void enable_depth_test(bool enabled);

///
/// Enable/disable culling
///
void enable_culling(bool enabled);

///
/// Enable/disable wireframe
///
void enable_wireframe(bool enabled);

///
/// Enable/disable edges. Attempts to avoid z-fighting but the method is not
/// robust.
///
void enable_edges(bool enabled, float tolerance = 0.9999f);

///
/// Enable/disable blending
///
void enable_blending(bool enabled);

///
/// Set blending to over operator
///
void set_blend_over();

///
/// Line width
///
void line_width(float w);

///
/// Set viewport
///
void set_viewport(const ym::vec4i& v);

// IMAGE FUNCTIONS -------------------------------------------------------------

///
/// Draw an texture tid of size img_w, img_h on a window of size win_w, win_h
/// with top-left corner at ox, oy with a zoom zoom.
///
void shade_image(
    uint tid, int win_w, int win_h, float ox, float oy, float zoom);

///
/// As above but includes an exposure/gamma correction.
///
void shade_image(uint tid, int win_w, int win_h, float ox, float oy, float zoom,
    ym::tonemap_type tmtype, float exposure, float gamma);

// TEXTURE FUNCTIONS -----------------------------------------------------------

///
/// Creates a texture with pixels values of size w, h with nc number of
/// components (1-4).
/// Internally use float if as_float and filtering if filter.
/// Returns the texture id.
///
uint make_texture(int w, int h, int nc, const float* pixels, bool linear,
    bool mipmap, bool as_float);

///
/// Creates a texture with pixels values of size w, h with nc number of
/// components (1-4).
/// Internally use srgb lookup if as_srgb and filtering if filter.
/// Returns the texture id.
///
uint make_texture(int w, int h, int nc, const unsigned char* pixels,
    bool linear, bool mipmap, bool as_srgb);

///
/// Updates the texture tid with new image data.
///
void update_texture(
    uint tid, int w, int h, int nc, const float* pixels, bool mipmap);

///
/// Updates the texture tid with new image data.
///
void update_texture(
    uint tid, int w, int h, int nc, const unsigned char* pixels, bool mipmap);

///
/// Creates a texture from an image.
/// Internally use float if as_float and filtering if filter.
/// Returns the texture id.
///
template <int N>
inline uint make_texture(const ym::image<ym::vec<float, N>>& img, bool linear,
    bool mipmap, bool as_float) {
    return make_texture(img.width(), img.height(), N, (float*)img.data(),
        linear, mipmap, as_float);
}

///
/// Creates a texture from an image.
/// Internally use srgb lookup if as_srgb and filtering if filter.
/// Returns the texture id.
///
template <int N>
inline uint make_texture(const ym::image<ym::vec<byte, N>>& img, bool linear,
    bool mipmap, bool as_srgb) {
    return make_texture(img.width(), img.height(), N, (byte*)img.data(), linear,
        mipmap, as_srgb);
}

///
/// Updates the texture tid with new image data.
///
template <int N>
inline void update_texture(
    uint tid, const ym::image<ym::vec<float, N>>& img, bool mipmap) {
    update_texture(
        tid, img.width(), img.height(), N, (float*)img.data(), mipmap);
}

///
/// Updates the texture tid with new image data.
///
template <int N>
inline void update_texture(
    uint tid, const ym::image<ym::vec<byte, N>>& img, bool mipmap) {
    update_texture(
        tid, img.width(), img.height(), N, (byte*)img.data(), mipmap);
}

///
/// Destroys the texture tid.
///
void clear_texture(uint* tid);

// BUFFER FUNCTIONS -----------------------------------------------------------

///
/// Creates a buffer with num elements of size size stored in values, where
/// content is dyanamic if dynamic.
/// Returns the buffer id.
///
uint make_buffer_raw(
    int size, const void* values, buffer_type btype, bool dynamic);

///
/// Creates a buffer with num elements of size size stored in values, where
/// content is dyanamic if dynamic.
/// Returns the buffer id.
///
template <typename T>
inline uint make_buffer(int num, const T* values, bool elements, bool dynamic) {
    return make_buffer_raw(num * sizeof(T), values,
        (elements) ? buffer_type::element : buffer_type::vertex, dynamic);
}

///
/// Creates a buffer with num elements of size size stored in values, where
/// content is dyanamic if dynamic.
/// Returns the buffer id.
///
template <typename T>
inline uint make_buffer(
    const std::vector<T>& values, bool elements, bool dynamic) {
    return make_buffer(values.size(), values.data(), elements, dynamic);
}

///
/// Creates a buffer with a block of data stored in values, where
/// content is dyanamic if dynamic.
/// Returns the buffer id.
///
template <typename T>
inline uint make_buffer(const T& values, buffer_type btype, bool dynamic) {
    return make_buffer_raw(sizeof(T), &values, btype, dynamic);
}

///
/// Updates the buffer bid with new data.
///
void update_buffer_raw(
    uint bid, int size, const void* values, buffer_type btype, bool dynamic);

///
/// Updates the buffer bid with new data.
///
template <typename T>
inline void update_buffer(
    uint bid, int num, const T* values, bool elements, bool dynamic) {
    update_buffer_raw(bid, num * sizeof(T), values,
        (elements) ? buffer_type::element : buffer_type::vertex, dynamic);
}

///
/// Updates the buffer bid with new data.
///
template <typename T>
inline void update_buffer(
    uint bid, const std::vector<T>& values, bool elements, bool dynamic) {
    update_buffer(bid, values.size(), values.data(), elements, dynamic);
}

///
/// Updates the buffer bid with new data.
///
template <typename T>
inline void update_buffer(
    uint bid, const T& values, buffer_type btype, bool dynamic) {
    update_buffer_raw(bid, sizeof(T), &values, btype, dynamic);
}

///
/// Destroys the buffer bid.
///
void clear_buffer(uint* bid);

// PROGRAM FUNCTIONS -----------------------------------------------------------

///
/// Creates and OpenGL vertex array object.
///
uint make_vertex_arrays();

///
/// Destroys the program pid and optionally the sahders vid and fid.
///
void clear_vertex_arrays(uint* aid);

///
/// Binds a vertex array object
///
void bind_vertex_arrays(uint aid);

///
/// Creates and OpenGL program from vertex and fragment code. Returns the
/// program id. Optionally return vertex and fragment shader ids. A VAO is
/// created.
///
uint make_program(const std::string& vertex, const std::string& fragment,
    uint* vid, uint* fid, uint* vao);

///
/// Destroys the program pid and optionally the sahders vid and fid.
///
void clear_program(uint* pid, uint* vid, uint* fid, uint* vao);

///
/// Binds a program
///
void bind_program(uint pid);

///
/// Get uniform location (simple GL wrapper that avoids GL includes)
///
int get_uniform_location(uint pid, const std::string& name);

///
/// Get uniform location (simple GL wrapper that avoids GL includes)
///
int get_uniform_block_index(uint pid, const std::string& name);

///
/// Get attrib location (simple GL wrapper that avoids GL includes)
///
int get_attrib_location(uint pid, const std::string& name);

///
/// Get the names of all uniforms
///
std::vector<std::pair<std::string, int>> get_uniforms_names(uint pid);

///
/// Get the names of all uniforms
///
inline std::unordered_map<std::string, int> get_uniforms_namemap(uint pid) {
    auto namemap = std::unordered_map<std::string, int>();
    for (auto kv : get_uniforms_names(pid)) namemap[kv.first] = kv.second;
    return namemap;
}

///
/// Get the names of all uniform blocks
///
std::vector<std::pair<std::string, int>> get_uniform_buffers_names(uint pid);

///
/// Get the names of all uniform blocks
///
inline std::unordered_map<std::string, int> get_uniform_buffers_namemap(
    uint pid) {
    auto namemap = std::unordered_map<std::string, int>();
    for (auto kv : get_uniform_buffers_names(pid))
        namemap[kv.first] = kv.second;
    return namemap;
}

///
/// Get the names of all attributes
///
std::vector<std::pair<std::string, int>> get_attributes_names(uint pid);

///
/// Get the names of all attributes
///
inline std::unordered_map<std::string, int> get_attributes_namemap(uint pid) {
    auto namemap = std::unordered_map<std::string, int>();
    for (auto kv : get_attributes_names(pid)) namemap[kv.first] = kv.second;
    return namemap;
}

///
/// Set uniform integer values val for program prog and variable loc.
/// The values have nc number of components (1-4) and count elements
/// (for arrays).
///
bool set_uniform_raw(uint prog, int var, const int* val, int nc, int count);

///
/// Set uniform float values val for program prog and variable var.
/// The values have nc number of components (1-4) and count elements
/// (for arrays).
///
bool set_uniform_raw(uint prog, int var, const float* val, int nc, int count);

///
/// Set uniform float values val for program prog and variable var.
///
inline bool set_uniform(uint prog, int var, bool val) {
    auto vali = (int)val;
    return set_uniform_raw(prog, var, &vali, 1, 1);
}

///
/// Set uniform float values val for program prog and variable var.
///
inline bool set_uniform(uint prog, int var, int val) {
    return set_uniform_raw(prog, var, &val, 1, 1);
}

///
/// Set uniform float values val for program prog and variable var.
///
inline bool set_uniform(uint prog, int var, float val) {
    return set_uniform_raw(prog, var, &val, 1, 1);
}

///
/// Set uniform float values val for program prog and variable var.
///
template <int N>
inline bool set_uniform(uint prog, int var, const ym::vec<float, N>& val) {
    return set_uniform_raw(prog, var, val.data(), N, 1);
}

///
/// Set uniform float values val for program prog and variable var.
///
template <int N>
inline bool set_uniform(uint prog, int var, const ym::vec<int, N>& val) {
    return set_uniform_raw(prog, var, val.data(), N, 1);
}

///
/// Set uniform float values val for program prog and variable var.
///
inline bool set_uniform(uint prog, int var, const ym::mat4f& val) {
    return set_uniform_raw(prog, var, (float*)val.data(), 16, 1);
}

///
/// Set uniform float values val for program prog and variable var.
///
inline bool set_uniform(uint prog, int var, const ym::frame3f& val) {
    return set_uniform_raw(prog, var, (float*)val.data(), 12, 1);
}

///
/// Set uniform float values val for program prog and variable var.
///
inline bool set_uniform(uint prog, int var, const int* val, int num) {
    return set_uniform_raw(prog, var, val, 1, num);
}

///
/// Set uniform float values val for program prog and variable var.
///
inline bool set_uniform(uint prog, int var, const float* val, int num) {
    return set_uniform_raw(prog, var, val, 1, num);
}

///
/// Set uniform float values val for program prog and variable var.
///
template <int N>
inline bool set_uniform(
    uint prog, int var, const ym::vec<float, N>* val, int num) {
    return set_uniform_raw(prog, var, val->data(), N, num);
}

///
/// Set uniform float values val for program prog and variable var.
///
template <int N>
inline bool set_uniform(
    uint prog, int var, const ym::vec<int, N>* val, int num) {
    return set_uniform_raw(prog, var, val->data(), N, num);
}

///
/// Set uniform float values val for program prog and variable var.
///
inline bool set_uniform(uint prog, int var, const ym::mat4f* val, int num) {
    return set_uniform_raw(prog, var, (float*)val, 16, num);
}

///
/// Set uniform float values val for program prog and variable var.
///
template <typename T>
inline bool set_uniform(uint prog, const std::string& var, const T& val) {
    auto loc = get_uniform_location(prog, var);
    if (loc < 0) return false;
    return set_uniform(prog, loc, val);
}

///
/// Set uniform float values val for program prog and variable var.
///
template <typename T>
inline bool set_uniform(
    uint prog, const std::string& var, const T* val, int num) {
    auto loc = get_uniform_location(prog, var);
    if (loc < 0) return false;
    return set_uniform(prog, loc, val, num);
}

///
/// Set uniform block bid for program prog and variable var.
///
bool set_uniform_block(uint prog, int var, int bid, uint bunit);

///
/// Set uniform texture id tid and unit tunit for program prog and variable var.
///
bool set_uniform_texture(
    uint prog, int var, const texture_info& tinfo, uint tunit);

///
/// Set uniform texture id tid and unit tunit for program prog and variable var.
/// Optionally sets the int variable varon to 0/1 whether the texture is enable
/// on not.
///
inline bool set_uniform_texture(
    uint prog, int var, int varon, const texture_info& tinfo, uint tunit) {
    if (!set_uniform_texture(prog, var, tinfo, tunit)) return false;
    if (!set_uniform(prog, varon, (bool)tinfo.txt_id)) return false;
    return true;
}

///
/// Set uniform texture id tid and unit tunit for program prog and variable var.
///
inline bool set_uniform_texture(
    uint prog, const std::string& var, const texture_info& tinfo, uint tunit) {
    auto loc = get_uniform_location(prog, var);
    if (loc < 0) return false;
    return set_uniform_texture(prog, loc, tinfo, tunit);
}

///
/// Set uniform texture id tid and unit tunit for program prog and variable var.
/// Optionally sets the int variable varon to 0/1 whether the texture is enable
/// on not.
///
inline bool set_uniform_texture(uint prog, const std::string& var,
    const std::string& varon, const texture_info& tinfo, uint tunit) {
    auto loc = get_uniform_location(prog, var);
    if (loc < 0) return false;
    auto locon = get_uniform_location(prog, varon);
    if (locon < 0) return false;
    return set_uniform_texture(prog, loc, locon, tinfo, tunit);
}

///
/// Sets a constant value for a vertex attribute for program prog and
/// variable var. The attribute has nc components.
///
bool set_vertattr_val_raw(uint prog, int var, const float* value, int nc);

///
/// Sets a constant value for a vertex attribute for program prog and
/// variable var. The attribute has nc components.
///
bool set_vertattr_val_raw(uint prog, int var, const int* value, int nc);

///
/// Sets a vartex attribute for program prog and variable var to the buffer bid.
/// The attribute has nc components and per-vertex values values.
///
bool set_vertattr_buffer(uint prog, int var, uint bid, int nc);

///
/// Sets a vartex attribute for program prog and variable var to the buffer bid.
/// The attribute has nc components and per-vertex values values.
///
bool set_vertattri_buffer(uint prog, int var, uint bid, int nc);

///
/// Sets a vartex attribute for program prog and variable var. The attribute
/// has nc components and either buffer bid or a single value def
/// (if bid is zero). Convenience wrapper to above functions.
///
bool set_vertattr_raw(uint prog, int var, uint bid, int nc, const float* def);

///
/// Sets a vartex attribute for program prog and variable var. The attribute
/// has nc components and either buffer bid or a single value def
/// (if bid is zero). Convenience wrapper to above functions.
///
bool set_vertattr_raw(uint prog, int var, uint bid, int nc, const int* def);

///
/// Sets a vartex attribute for program prog and variable var. The attribute
/// is either a buffer bid or a single value def
/// (if bid is zero). Convenience wrapper to above functions.
///
template <int N>
inline bool set_vertattr(
    uint prog, int var, uint bid, const ym::vec<float, N>& def) {
    return set_vertattr_raw(prog, var, bid, N, def.data());
}

///
/// Sets a vartex attribute for program prog and variable var. The attribute
/// is either a buffer bid or a single value def
/// (if bid is zero). Convenience wrapper to above functions.
///
template <int N>
inline bool set_vertattr(
    uint prog, int var, uint bid, const ym::vec<int, N>& def) {
    return set_vertattr_raw(prog, var, bid, N, def.data());
}

///
/// Sets a vartex attribute for program prog and variable var. The attribute
/// is either a buffer bid or a single value def
/// (if bid is zero). Convenience wrapper to above functions.
///
template <int N>
inline bool set_vertattr(
    uint prog, const std::string& var, uint bid, const ym::vec<float, N>& def) {
    auto loc = get_attrib_location(prog, var);
    if (loc < 0) return false;
    return set_vertattr(prog, loc, bid, def);
}

///
/// Sets a vartex attribute for program prog and variable var. The attribute
/// is either a buffer bid or a single value def
/// (if bid is zero). Convenience wrapper to above functions.
///
template <int N>
inline bool set_vertattr(
    uint prog, const std::string& var, uint bid, const ym::vec<int, N>& def) {
    auto loc = get_attrib_location(prog, var);
    if (loc < 0) return false;
    return set_vertattr(prog, loc, bid, def);
}

///
/// Draws nelems elements elem of type etype.
///
bool draw_elems(int nelems, uint bid, etype etype);

// STANDARD SHADER FUNCTIONS ---------------------------------------------------

///
/// Shade with a physically-based standard shader based on Phong/GGX.
///
namespace stdshader {
///
/// Initialize a standard shader.
///
uint make_program(uint* vao);

///
/// Starts a frame by setting exposure/gamma values, camera transforms and
/// projection. Sets also whether to use full shading or a quick eyelight
/// preview.
///
void begin_frame(uint prog, uint vao, bool shade_eyelight, float img_exposure,
    ym::tonemap_type img_tonemap, float img_gamma,
    const ym::mat4f& camera_xform, const ym::mat4f& camera_xform_inv,
    const ym::mat4f& camera_proj);

///
/// Ends a frame.
///
void end_frame();

///
/// Set num lights with position pos, color ke, type ltype. Also set the
/// ambient illumination amb.
///
void set_lights(uint prog, const ym::vec3f& amb, int num, ym::vec3f* pos,
    ym::vec3f* ke, ltype* ltype);

///
/// Begins drawing a shape with transform xform.
///
void begin_shape(uint prog, const ym::mat4f& xform);

///
/// End shade drawing.
///
void end_shape();

///
/// Set material values for emission only (constant color).
/// Indicates textures ids with the correspoinding XXX_txt variables.
/// Works for points/lines/triangles. Element type set by draw_XXX calls.
///
void set_material_emission_only(uint prog, etype etype, const ym::vec3f& ke,
    float op, const texture_info& ke_txt, bool double_sided, bool alpha_cutout);

///
/// Set material values with emission ke, diffuse kd, specular ks and
/// specular roughness rs, opacity op. Indicates textures ids with the
/// correspoinding XXX_txt variables. Sets also normal and occlusion
/// maps. Works for points/lines/triangles (diffuse for points,
/// Kajiya-Kay for lines, GGX/Phong for triangles). Element type set by
/// draw_XXX calls.
///
void set_material_generic(uint prog, etype etype, const ym::vec3f& ke,
    const ym::vec3f& kd, const ym::vec3f& ks, float rs, float op,
    const texture_info& ke_txt, const texture_info& kd_txt,
    const texture_info& ks_txt, const texture_info& rs_txt,
    const texture_info& norm_txt, const texture_info& occ_txt, bool use_phong,
    bool double_sided, bool alpha_cutout);

///
/// Set material values for glTF specular-roughness PBR shader,
/// with emission ke, base color kb, opacity op, metallicity km and
/// specular roughness rs. Uses basecolor-opacity texture kb_txt and
/// metallic-roughness texture km_txt. Sets also normal and occlusion
/// maps. Works for points/lines/triangles (diffuse for points, Kajiya-Kay
/// for lines, GGX/Phong for triangles). Element type set by draw_XXX calls.
///
void set_material_gltf_metallic_roughness(uint prog, etype etype,
    const ym::vec3f& ke, const ym::vec3f& kb, float km, float rs, float op,
    const texture_info& ke_txt, const texture_info& kb_txt,
    const texture_info& km_txt, const texture_info& norm_txt,
    const texture_info& occ_txt, bool use_phong, bool double_sided,
    bool alpha_cutout);

///
/// Set material values for glTF specular-roughness PBR shader,
/// with emission ke, diffuse color kd, opacity op, specular ks and
/// specular glossiness rs. Uses diffuse-opacity texture kd_txt and
/// specular-glpossiness texture ks_txt. Sets also normal and occlusion
/// maps. Works for points/lines/triangles (diffuse for points, Kajiya-Kay
/// for lines, GGX/Phong for triangles). Element type set by draw_XXX calls.
///
void set_material_gltf_specular_glossiness(uint prog, etype etype,
    const ym::vec3f& ke, const ym::vec3f& kd, const ym::vec3f& ks, float rs,
    float op, const texture_info& ke_txt, const texture_info& kd_txt,
    const texture_info& ks_txt, const texture_info& norm_txt,
    const texture_info& occ_txt, bool use_phong, bool double_sided,
    bool alpha_cutout);

///
/// Convertes a phong exponent to roughness.
///
float specular_exponent_to_roughness(float n);

///
/// Set vertex data with position pos, normals norm, texture coordinates
/// texcoord and per-vertex color color and tangent space tangsp.
///
void set_vert(uint prog, const ym::vec3f* pos, const ym::vec3f* norm,
    const ym::vec2f* texcoord, const ym::vec4f* color, const ym::vec4f* tangsp);

///
/// Set vertex data with buffers for position pos, normals norm, texture
/// coordinates texcoord, per-vertex color color and tangent space tangsp.
///
void set_vert(
    uint prog, uint pos, uint norm, uint texcoord, uint color, uint tangsp);

///
/// Set vertex data with buffers for skinning.
///
void set_vert_skinning(
    uint prog, uint weights, uint joints, int nxforms, const ym::mat4f* xforms);

///
/// Set vertex data with buffers for skinning.
///
void set_vert_gltf_skinning(
    uint prog, uint weights, uint joints, int nxforms, const ym::mat4f* xforms);

///
/// Disables vertex skinning.
///
void set_vert_skinning_off(uint prog);

///
/// Set the object as highlight color.
///
void set_highlight(uint prog, const ym::vec4f& highlight);

///
/// Draw num elements elem of type etype.
///
inline void draw_elems(uint prog, int num, uint bid, etype etype) {
    yglu::draw_elems(num, bid, etype);
}

}  // namespace stdshader

}  // namespace yglu

#endif
