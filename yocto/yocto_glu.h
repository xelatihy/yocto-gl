///
/// YOCTO_GLU: a set of utilities to draw on screen with OpenGL 2.1. Mostly
/// used to make quick viewers. Not sure it is helpful to others (OpenGL is
/// really a portability mess, but there is nothing else). The library is split
/// into two namespaces, legacy and modern, to indicate functions for fixed
/// function pipeline and for modern shaders.
///
///
/// FEATURES:
///
/// 1. different OpenGL versions are supported in the legacy and modern
/// namespaces, resectively 1.X and 3.X
/// 2. image viewing with legacy::draw_image() or modern::shade_image()
/// - the modern interface supports exposure/gamma correction
/// 3. texture utilies to quickly create/update textures
/// - create textures with legacy::make_texture() or modern::make_texture()
/// - create textures with legacy::make_texture() or modern::update_texture()
/// - delete textures with legacy::clear_texture() or monder::clear_texture()
/// 4. program utilities in modern namespace
/// - buffer objects with create_buffer(), update_buffer()
/// - program creation/cleaning with make_program(), clear_program()
/// - uniforms with set_uniform()
/// - vertex attrib with set_vertattr() and set_vertattr_buff()
/// - draw_elems()
/// 5. a standard shader for GGX fragment shading and multiple lights in
/// the stdshader namespace
/// - initialize the shaders with make_program()
/// - start/end each frame with begin_frame(), end_frame()
/// - define lights with set_lights(light params)
/// - start/end each shape with begin_shape(), end_shape()
/// - define material parameters with set_material()
/// - define vertices with set_vert()
/// - draw elements with draw_elems(element data)
/// 6. standard drawing methods legacy namespace with interfaces as in 6.
/// 7. user interface functions based on GLFW/Nuklear
///
/// The interface for each function is described in details in the interface
/// section of this file.
///
/// COMPILATION:
///
/// To use the library include the .h and compile the .cpp. To use this library
/// as a header-only library, define YBVH_INLINE before including this file.
/// The library depends on GLEW for OpenGL functions, GLFW and Nuklear for UI
/// function. The latter needed to be avilable locally. To disable compilation
/// of UI function use YGLU_NO_GLFW and YGLU_NO_NUKLEAR.
///
///
/// HISTORY:
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

// compilation options
#ifdef YGLU_INLINE
#define YGLU_API inline
#else
#define YGLU_API
#endif

#include <array>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// INTERFACE
// -----------------------------------------------------------------------------

namespace yglu {

//
// Typedefs for vec/mat types
//
using float2 = std::array<float, 2>;
using float3 = std::array<float, 3>;
using float4 = std::array<float, 4>;
using float4x4 = std::array<std::array<float, 4>, 4>;
using int2 = std::array<int, 2>;
using int3 = std::array<int, 3>;
using int4 = std::array<int, 4>;
using byte4 = std::array<unsigned char, 4>;

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
    quad = 4
};

///
/// Light types
///
enum struct ltype : int {
    point = 0,       ///< point lights
    directional = 1  ///< directional lights
};

///
/// Shortcut for GLuint.
///
/// Would be better to use GLuint but this avoids issues with including GL
/// headers.
///
using uint = unsigned int;

///
/// Checks for GL error and then prints
///
YGLU_API bool check_error(bool print = true);

///
/// Clear window
///
YGLU_API void clear_buffers(const float4& background = {0, 0, 0, 0});

///
/// Enable/disable depth test
///
YGLU_API void enable_depth_test(bool enabled);

///
/// Enable/disable culling
///
YGLU_API void enable_culling(bool enabled);

///
/// Enable/disable wireframe
///
YGLU_API void enable_wireframe(bool enabled);

///
/// Enable/disable edges. Attempts to avoid z-fighting but the method is not
/// robust.
///
YGLU_API void enable_edges(bool enabled, float tolerance = 0.9999f);

///
/// Line width
///
YGLU_API void line_width(float w);

///
/// Set viewport
///
YGLU_API void set_viewport(const int4& v);

// -----------------------------------------------------------------------------
// LEGACY FUNCTIONS
// -----------------------------------------------------------------------------

///
/// Legacy Functions (OpenGL 1.X)
///
namespace legacy {

// IMAGE FUNCTIONS -------------------------------------------------------------

///
/// Draw an texture tid of size img_w, img_h on a window of size win_w, win_h
/// with top-left corner at ox, oy with a zoom zoom.
///
YGLU_API void draw_image(uint tid, int img_w, int img_h, int win_w, int win_h,
    float ox, float oy, float zoom);

///
/// Reads an image back to memory.
///
YGLU_API void read_imagef(float* pixels, int w, int h, int nc);

///
/// Creates a texture with pixels values of size w, h with nc number of
/// components (1-4).
/// Internally use filtering if filter.
/// Returns the texture id.
///
YGLU_API uint make_texture(
    int w, int h, int nc, const float* pixels, bool linear, bool mipmap);

///
/// Creates a texture with pixels values of size w, h with nc number of
/// components (1-4).
/// Internally use filtering if filter.
/// Returns the texture id.
///
YGLU_API uint make_texture(int w, int h, int nc, const unsigned char* pixels,
    bool linear, bool mipmap);

///
/// Updates the texture tid with new image data.
///
YGLU_API void update_texture(
    uint tid, int w, int h, int nc, const float* pixels, bool mipmap);

///
/// Updates the texture tid with new image data.
///
YGLU_API void update_texture(
    uint tid, int w, int h, int nc, const unsigned char* pixels, bool mipmap);

///
/// Destroys the texture tid.
///
YGLU_API void clear_texture(uint* tid);

///
/// Starts a frame by setting exposure/gamma values, camera transforms and
/// projection. Sets also whether to use full shading or a quick eyelight
/// preview.
/// If eyelight is disabled, sets num lights with position pos,
/// color ke, type ltype. Also set the ambient illumination amb.
///
YGLU_API void begin_frame(const float4x4& camera_xform,
    const float4x4& camera_xform_inv, const float4x4& camera_proj,
    bool eyelight, bool scale_kx);

///
/// Ends a frame.
///
YGLU_API void end_frame();

///
/// Set num lights with position pos, color ke, type ltype. Also set the ambient
/// illumination amb.
///
YGLU_API void set_lights(const float3& amb, int num, const float3* pos,
    const float3* ke, const ltype* ltype);

///
/// Begins drawing a shape with transform xform.
///
YGLU_API void begin_shape(const float4x4& xform);

///
/// End shade drawing.
///
YGLU_API void end_shape();

///
/// Set material values with emission ke, diffuse kd, specular ks and
/// specular exponent ns. Indicates textures ids with the correspoinding XXX_txt
/// variables.
///
YGLU_API void set_material(const float3& ke, const float3& kd, const float3& ks,
    float ns, int kd_txt, bool scale_kx);

///
/// Convertes a phong exponent to roughness.
///
YGLU_API float specular_roughness_to_exponent(float r);

///
/// Draw num elements elem of type etype.
///
YGLU_API void draw_elems(int num, const int* elem, etype etype,
    const float3* pos, const float3* norm, const float2* texcoord,
    const float3* color);

///
/// Draw num elements elem of type etype.
///
YGLU_API void draw_points(int num, const int* elem, const float3* pos,
    const float3* norm, const float2* texcoord, const float3* color);

///
/// Draw num elements elem of type etype.
///
YGLU_API void draw_lines(int num, const int2* elem, const float3* pos,
    const float3* norm, const float2* texcoord, const float3* color);

///
/// Draw num elements elem of type etype.
///
YGLU_API void draw_triangles(int num, const int3* elem, const float3* pos,
    const float3* norm, const float2* texcoord, const float3* color);

}  // namespace

// -----------------------------------------------------------------------------
// MODERN FUNCTIONS
// -----------------------------------------------------------------------------

///
/// Modern OpenGL (OpenGL > 3.2)
///
namespace modern {

// IMAGE FUNCTIONS -------------------------------------------------------------

///
/// Draw an texture tid of size img_w, img_h on a window of size win_w, win_h
/// with top-left corner at ox, oy with a zoom zoom.
///
YGLU_API void shade_image(uint tid, int img_w, int img_h, int win_w, int win_h,
    float ox, float oy, float zoom);

///
/// As above but includes an exposure/gamma correction.
///
YGLU_API void shade_image(uint tid, int img_w, int img_h, int win_w, int win_h,
    float ox, float oy, float zoom, float exposure, float gamma_);

// TEXTURE FUNCTIONS -----------------------------------------------------------

///
/// Creates a texture with pixels values of size w, h with nc number of
/// components (1-4).
/// Internally use float if as_float and filtering if filter.
/// Returns the texture id.
///
YGLU_API uint make_texture(int w, int h, int nc, const float* pixels,
    bool linear, bool mipmap, bool as_float);

///
/// Creates a texture with pixels values of size w, h with nc number of
/// components (1-4).
/// Internally use srgb lookup if as_srgb and filtering if filter.
/// Returns the texture id.
///
YGLU_API uint make_texture(int w, int h, int nc, const unsigned char* pixels,
    bool linear, bool mipmap, bool as_srgb);

///
/// Updates the texture tid with new image data.
///
YGLU_API void update_texture(
    uint tid, int w, int h, int nc, const float* pixels, bool mipmap);

///
/// Updates the texture tid with new image data.
///
YGLU_API void update_texture(
    uint tid, int w, int h, int nc, const unsigned char* pixels, bool mipmap);

///
/// Destroys the texture tid.
///
YGLU_API void clear_texture(uint* tid);

// BUFFER FUNCTIONS -----------------------------------------------------------

///
/// Creates a buffer with num elements of size size stored in values, where
/// content is dyanamic if dynamic.
/// Returns the buffer id.
///
YGLU_API uint make_buffer(
    int num, int size, const void* values, bool elements, bool dynamic);

///
/// Updates the buffer bid with new data.
///
YGLU_API void update_buffer(uint bid, int num, int size, const void* values,
    bool elements, bool dynamic);

///
/// Destroys the buffer bid.
///
YGLU_API void clear_buffer(uint* bid);

// PROGRAM FUNCTIONS -----------------------------------------------------------

///
/// Creates and OpenGL vertex array object.
///
YGLU_API uint make_vertex_arrays();

///
/// Destroys the program pid and optionally the sahders vid and fid.
///
YGLU_API void clear_vertex_arrays(uint* aid);

///
/// Creates and OpenGL program from vertex and fragment code. Returns the
/// program id. Optionally return vertex and fragment shader ids. A VAO has to
/// be
/// bound before this.
///
YGLU_API uint make_program(const std::string& vertex,
    const std::string& fragment, uint* vid, uint* fid);

///
/// Destroys the program pid and optionally the sahders vid and fid.
///
YGLU_API void clear_program(uint* pid, uint* vid, uint* fid);

///
/// Set uniform integer values val for program prog and variable var.
/// The values have nc number of components (1-4) and count elements
/// (for arrays).
///
YGLU_API bool set_uniform(
    uint prog, const std::string& var, const int* val, int nc, int count);

///
/// Set uniform float values val for program prog and variable var.
/// The values have nc number of components (1-4) and count elements
/// (for arrays).
///
YGLU_API bool set_uniform(
    uint prog, const std::string& var, const float* val, int nc, int count);

///
/// Set uniform texture id tid and unit tunit for program prog and variable var.
/// Optionally sets the int variable varon to 0/1 whether the texture is enable
/// on not.
///
YGLU_API bool set_uniform_texture(uint prog, const std::string& var,
    const std::string& varon, uint tid, uint tunit);

///
/// Sets a constant value for a vertex attribute for program prog and
/// variable var. The attribute has nc components.
///
YGLU_API bool set_vertattr_val(
    uint prog, const std::string& var, const float* value, int nc);

///
/// Sets a vartex attribute for program prog and variable var to the buffer bid.
/// The attribute has nc components and per-vertex values values.
///
YGLU_API bool set_vertattr_buffer(
    uint prog, const std::string& var, uint bid, int nc);

///
/// Sets a vartex attribute for program prog and variable var. The attribute
/// has nc components and either buffer bid or a single value def
/// (if bid is zero). Convenience wrapper to above functions.
///
YGLU_API bool set_vertattr(
    uint prog, const std::string& var, uint bid, int nc, const float* def);

///
/// Draws nelems elements elem of type etype.
///
YGLU_API bool draw_elems(int nelems, uint bid, etype etype);

}  // namespace

// STANDARD SHADER FUNCTIONS ---------------------------------------------------

///
/// Shade with a physically-based standard shader based on Phong/GGX.
///
namespace stdshader {

///
/// Initialize a standard shader.
///
YGLU_API void make_program(uint* pid, uint* vao);

///
/// Tone mapping presets
///
enum struct tonemap_type {
    /// default
    def = 0,
    /// linear
    linear = 1,
    /// srgb
    srgb = 2,
    /// linear
    gamma = 3,
    /// filmic
    filmic = 4
};

///
/// Starts a frame by setting exposure/gamma values, camera transforms and
/// projection. Sets also whether to use full shading or a quick eyelight
/// preview.
///
YGLU_API void begin_frame(uint prog, uint vao, bool shade_eyelight,
    float img_exposure, tonemap_type img_tonemap, float img_gamma,
    const float4x4& camera_xform, const float4x4& camera_xform_inv,
    const float4x4& camera_proj);

///
/// Ends a frame.
///
YGLU_API void end_frame();

///
/// Set num lights with position pos, color ke, type ltype. Also set the ambient
/// illumination amb.
///
YGLU_API void set_lights(uint prog, const float3& amb, int num, float3* pos,
    float3* ke, ltype* ltype);

///
/// Begins drawing a shape with transform xform.
///
YGLU_API void begin_shape(uint prog, const float4x4& xform);

///
/// End shade drawing.
///
YGLU_API void end_shape();

///
/// Set material values with emission ke, diffuse kd, specular ks and
/// specular roughness rs. Indicates textures ids with the correspoinding
/// XXX_txt
/// variables. Uses GGX by default, but can switch to Phong is needed. Works
/// for points, lines and triangle/quads based on etype.
///
YGLU_API void set_material(uint prog, const float3& ke, const float3& kd,
    const float3& ks, float rs, int ke_txt, int kd_txt, int ks_txt, int rs_txt,
    bool use_phong);

///
/// Convertes a phong exponent to roughness.
///
YGLU_API float specular_exponent_to_roughness(float n);

///
/// Set vertex data with position pos, normals norm, texture coordinates
/// texcoord
/// and per-vertex color color.
///
YGLU_API void set_vert(uint prog, const float3* pos, const float3* norm,
    const float2* texcoord, const float3* color);

///
/// Set vertex data with buffers for position pos, normals norm, texture
/// coordinates texcoord and per-vertex color color.
///
YGLU_API void set_vert(
    uint prog, uint pos, uint norm, uint texcoord, uint color);

///
/// Draw num elements elem of type etype.
///
YGLU_API void draw_elems(uint prog, int num, uint bid, etype etype);

///
/// Draw num elements elem of type etype.
///
YGLU_API void draw_points(uint prog, int num, uint bid);

///
/// Draw num elements elem of type etype.
///
YGLU_API void draw_lines(uint prog, int num, uint bid);

///
/// Draw num elements elem of type etype.
///
YGLU_API void draw_triangles(uint prog, int num, uint bid);

}  // namespace

// -----------------------------------------------------------------------------
// USER INTERFACE FUNCTIONS
// -----------------------------------------------------------------------------

#ifndef YGLU_NO_GLFW

///
/// Modern OpenGL (OpenGL > 3.2)
///
namespace ui {
///
/// Forward declaration of window
///
struct window;

///
/// Key callback
///
typedef void (*text_callback)(window*, unsigned int);

///
/// Window refresh callback
///
typedef void (*refresh_callback)(window*);

///
/// initialize glfw
///
window* init_window(int width, int height, const std::string& title,
    bool legacy_gl, void* user_pointer = nullptr);

///
/// Clear glfw
///
void clear_window(window* window);

///
/// initialize glfw
///
void set_callbacks(
    window* win, text_callback text_cb, refresh_callback refresh_cb = nullptr);

///
/// Set window title
///
void set_window_title(window* win, const std::string& title);

///
/// Wait events
///
void wait_events(window* win);

///
/// Poll events
///
void poll_events(window* win);

///
/// Swap buffers
///
void swap_buffers(window* win);

///
/// Should close
///
bool should_close(window* win);

///
/// User pointer
///
void* get_user_pointer(window* window);

///
/// Mouse button
///
int get_mouse_button(window* window);

///
/// Mouse position
///
int2 get_mouse_posi(window* window);

///
/// Mouse position
///
float2 get_mouse_posf(window* window);

///
/// Window size
///
int2 get_window_size(window* window);

///
/// Framebuffer size
///
int2 get_framebuffer_size(window* window);

///
/// Read pixels
///
std::vector<byte4> get_screenshot(
    window* win, int2& wh, bool flipy = true, bool back = false);

#ifndef YGLU_NO_NUKLEAR

///
/// Init ui
///
void init_widgets(window* win);

///
/// Clear ui
///
void clear_widgets(window* win);

///
/// Begin draw widget
///
bool begin_widgets(window* win);

///
/// End draw widget
///
void end_widgets(window* win);

///
/// Dynamic layout for next widgets
///
void dynamic_widget_layout(window* win, int n);

///
/// Label widget
///
void label_widget(window* win, const std::string& lbl);

///
/// Label and int widget
///
void int_label_widget(window* win, const std::string& lbl, int val);

///
/// Label and float widget
///
void float_label_widget(window* win, const std::string& lbl, float val);

///
/// Int widget
///
void int_widget(window* win, const std::string& lbl, int* val, int min, int max,
    int incr = 1);

///
/// Float widget
///
void float_widget(window* win, const std::string& lbl, float* val, float min,
    float max, float incr = 1.0f);

///
/// Bool widget
///
void bool_widget(window* win, const std::string& lbl, bool* val);

///
/// Enum widget
///
void enum_widget(window* win, const std::string& lbl, int* val,
    const std::vector<std::pair<std::string, int>>& labels);

///
/// Button widget
///
bool button_widget(window* win, const std::string& lbl);

///
/// Whether widget are active
///
bool get_widget_active(window* win);

#endif
}  // namespace

#endif

}  // namespace

// -----------------------------------------------------------------------------
// INCLUDE FOR HEADER-ONLY MODE
// -----------------------------------------------------------------------------

#ifdef YGLU_INLINE
#include "yocto_glu.cpp"
#endif

#endif
