//
// # Yocto/GLUtils: Tiny OpenGL utilities for writing simple interactive apps.
//
// Small set of utilities to support writing OpenGL 3.3, manage
// windows with GLFW and draw immediate-mode widgets with ImGui.
//
// 1. texture and buffer objects with `gltexture` and `glbuffer`
//     - create textures/buffers with appropriate constructors
//     - check validity with `is_valid()`
//     - update textures/buffers with `update()` functions
//     - delete textures/buffers with `clear()`
//     - bind/unbind textures/buffers with `bind()`/`unbind()`
//     - draw elements with `gl_buffer::draw_elems()`
// 2. program objects with `glprogram`
//     - program creation with constructor
//     - check validity with `is_valid()`
//     - delete with `clear()`
//     - uniforms with `set_gluniform()`
//     - vertex attrib with `set_glattribute()`
//     - draw elements with `gl_buffer::draw_elems()`
// 3. image viewing with `glimage_program`, with support for tone mapping.
// 4. draw surfaces and hair with GGX/Kayjia-Kay with `glsurface_program`
//     - initialize the program with constructor
//     - check validity with `is_valid()`
//     - start/end each frame with `begin_frame()`, `end_frame()`
//     - define lights with `set_lights()`
//     - start/end each shape with `begin_shape()`, `end_shape()`
//     - define material Parameters with `set_material()`
//     - define vertices with `set_vert()`
//     - draw elements with `draw_elems()`
// 5. also includes other utlities for quick OpenGL hacking
// 6. GLFW window with `glwindow`
//     - create with constructor
//     - delete with `clear()`
//     - set callbacks with `set_callbacks()`
//     - includes carious utilities to query window, mouse and keyboard
// 7. immediate mode widgets using ImGui
//     - init with `init_widget()`
//     - use the various widget calls to draw the widget and handle events
//
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2018 Fabio Pellacini
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
//

#ifndef _YGLU_H_
#define _YGLU_H_

// -----------------------------------------------------------------------------
// COMPILATION OPTIONS AND INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_math.h"

#include <functional>
#include <map>

// -----------------------------------------------------------------------------
// OPENGL OBJECTS AND FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// OpenGL shape element types.
enum struct glelem_type : int { point = 1, line = 2, triangle = 3 };

// OpenGL light types.
enum struct gllight_type : int {
    point = 0,
    directional = 1,
};

// OpenGL lights
struct gllights {
    std::vector<vec3f> pos;          // light positions
    std::vector<vec3f> ke;           // light intensities
    std::vector<gllight_type> type;  // light types
};

// Checks for GL error and then prints.
bool check_glerror(bool print = true);

// Clear window.
void clear_glbuffers(const vec4f& background = {0, 0, 0, 0});

// Enable/disable depth test, culling, wireframe and blending.
void enable_gldepth_test(bool enabled);
void enable_glculling(bool enabled, bool front = false, bool back = true);
void enable_glwireframe(bool enabled);
void enable_glblending(bool enabled);
void set_glblend_over();

// Set viewport.
void set_glviewport(const vec4i& v);
void set_glviewport(const vec2i& v);

// Reads an image from the the framebuffer.
void read_glimagef(float* pixels, int w, int h, int nc);

// OpenGL texture object. Members are not part of the public API.
struct gltexture {
    uint tid = 0;           // texture id
    int width = 0;          // width
    int height = 0;         // height
    int ncomp = 0;          // number of components
    bool as_float = false;  // stored as float
    bool as_srgb = true;    // stored as sRGB
    bool mipmap = true;     // store with mipmaps
    bool linear = true;     // use linear interpolation
};

// Implementation of update_texture.
void update_gltexture(gltexture& txt, int w, int h, int nc, const void* pixels,
    bool floats, bool linear, bool mipmap, bool as_float, bool as_srgb);

// Updates a texture with pixels values of size w, h with nc number of
// components (1-4). Internally use bytes/floats (as_float), linear/sRGB
// (as_srgb) nearest/linear filtering (linear) and mipmmapping (mipmap).
inline void update_gltexture(gltexture& txt, int w, int h, int nc,
    const float* pixels, bool linear, bool mipmap, bool as_float) {
    update_gltexture(
        txt, w, h, nc, pixels, true, linear, mipmap, as_float, false);
}
inline void update_gltexture(gltexture& txt, int w, int h, int nc,
    const unsigned char* pixels, bool linear, bool mipmap, bool as_srgb) {
    update_gltexture(
        txt, w, h, nc, pixels, false, linear, mipmap, false, as_srgb);
}

// Updates a texture with pixels values from an image.
inline void update_gltexture(gltexture& txt, int width, int height,
    const std::vector<vec4b>& pixels, bool linear, bool mipmap, bool as_srgb) {
    update_gltexture(txt, width, height, 4, (unsigned char*)pixels.data(),
        linear, mipmap, as_srgb);
}
inline void update_gltexture(gltexture& txt, int width, int height,
    const std::vector<vec4f>& pixels, bool linear, bool mipmap, bool as_float) {
    update_gltexture(
        txt, width, height, 4, (float*)pixels.data(), linear, mipmap, as_float);
}

// Binds/unbinds a texture to a texture unit.
void bind_gltexture(const gltexture& txt, uint unit);
void unbind_gltexture(const gltexture& txt, uint unit);

// Clears the texture.
void clear_gltexture(gltexture& txt);

// Get texture id and check if defined.
inline uint get_gltexture_id(const gltexture& txt) { return txt.tid; }
inline bool is_gltexture_valid(const gltexture& txt) { return (bool)txt.tid; }

// Wrap values for OpenGL texture.
enum struct gltexture_wrap { not_set, repeat, clamp, mirror };

// Filter values for OpenGL texture.
enum struct gltexture_filter {
    not_set,
    linear,
    nearest,
    linear_mipmap_linear,
    nearest_mipmap_nearest,
    linear_mipmap_nearest,
    nearest_mipmap_linear
};

// OpenGL texture parameters.
struct gltexture_info {
    gltexture txt = {};                               // texture
    int texcoord = 0;                                 // texture coordinate set
    float scale = 1;                                  // texture scale
    gltexture_wrap wrap_s = gltexture_wrap::not_set;  // wrap mode s
    gltexture_wrap wrap_t = gltexture_wrap::not_set;  // wrap mode s
    gltexture_filter filter_mag = gltexture_filter::not_set;  // mag filter
    gltexture_filter filter_min = gltexture_filter::not_set;  // min filter
};

// OpenGL vertex/element buffer. Members are not part of the public API.
struct glvertex_buffer {
    uint bid = 0;        // buffer id
    int num = 0;         // number of elements
    int ncomp = 0;       // number of components
    bool elems = false;  // element buffer
};

// Updates vertex/element buffers of floats/ints respectively.
void update_glbuffer(glvertex_buffer& buf, bool elems, int num, int ncomp,
    const void* values, bool dynamic);

// Updates the buffer with new data.
inline void update_glbuffer(glvertex_buffer& buf, bool elems,
    const std::vector<float>& values, bool dynamic = false) {
    update_glbuffer(buf, elems, values.size(), 1, values.data(), dynamic);
}
inline void update_glbuffer(glvertex_buffer& buf, bool elems,
    const std::vector<vec2f>& values, bool dynamic = false) {
    update_glbuffer(
        buf, elems, values.size(), 2, (const float*)values.data(), dynamic);
}
inline void update_glbuffer(glvertex_buffer& buf, bool elems,
    const std::vector<vec3f>& values, bool dynamic = false) {
    update_glbuffer(
        buf, elems, values.size(), 3, (const float*)values.data(), dynamic);
}
inline void update_glbuffer(glvertex_buffer& buf, bool elems,
    const std::vector<vec4f>& values, bool dynamic = false) {
    update_glbuffer(
        buf, elems, values.size(), 4, (const float*)values.data(), dynamic);
}
inline void update_glbuffer(glvertex_buffer& buf, bool elems,
    const std::vector<int>& values, bool dynamic = false) {
    update_glbuffer(buf, elems, values.size(), 1, values.data(), dynamic);
}
inline void update_glbuffer(glvertex_buffer& buf, bool elems,
    const std::vector<vec2i>& values, bool dynamic = false) {
    update_glbuffer(
        buf, elems, values.size(), 2, (const int*)values.data(), dynamic);
}
inline void update_glbuffer(glvertex_buffer& buf, bool elems,
    const std::vector<vec3i>& values, bool dynamic = false) {
    update_glbuffer(
        buf, elems, values.size(), 3, (const int*)values.data(), dynamic);
}
inline void update_glbuffer(glvertex_buffer& buf, bool elems,
    const std::vector<vec4i>& values, bool dynamic = false) {
    update_glbuffer(
        buf, elems, values.size(), 4, (const int*)values.data(), dynamic);
}

// Binds/unbinds the buffer at a particular attribute location.
void bind_glbuffer(const glvertex_buffer& buf, uint vattr);
void unbind_glbuffer(const glvertex_buffer& buf, uint vattr);
void unbind_glbuffer(uint vattr);

// Get buffer id and if valid.
inline uint get_glbuffer_id(const glvertex_buffer& buf) { return buf.bid; }
inline bool is_glbuffer_valid(const glvertex_buffer& buf) {
    return (bool)buf.bid;
}
inline bool is_glbuffer_empty(const glvertex_buffer& buf) {
    return !buf.bid || !buf.num;
}

// Clears OpenGL state.
void clear_glbuffer(glvertex_buffer& buf);

// Draws elements.
void draw_glelems(const glvertex_buffer& buf);

// OpenGL program. Members are not part of the public API.
struct glprogram {
    uint pid = 0;  // program id
    uint vid = 0;  // vertex shader is
    uint fid = 0;  // fragment shader id
    uint vao = 0;  // vertex array object id
};

// Creates an OpenGL program from vertex and fragment code.
glprogram make_glprogram(
    const std::string& vertex, const std::string& fragment);

// Get uniform and attribute locations.
int get_gluniform_location(const glprogram& prog, const std::string& name);
int get_glattrib_location(const glprogram& prog, const std::string& name);

// Set uniform values.
void set_gluniform(const glprogram& prog, int var, bool val);
void set_gluniform(const glprogram& prog, int var, int val);
void set_gluniform(const glprogram& prog, int var, float val);
void set_gluniform(const glprogram& prog, int var, const vec2f& val);
void set_gluniform(const glprogram& prog, int var, const vec3f& val);
void set_gluniform(const glprogram& prog, int var, const vec4f& val);
void set_gluniform(const glprogram& prog, int var, const mat4f& val);
void set_gluniform(const glprogram& prog, int var, const frame3f& val);
template <typename T>
inline void set_gluniform(
    const glprogram& prog, const std::string& var, const T& val) {
    set_gluniform(prog, get_gluniform_location(prog, var), val);
}

// Set uniform texture.
void set_gluniform_texture(
    const glprogram& prog, int pos, const gltexture_info& tinfo, uint tunit);
// Set uniform texture with an additionasl texture enable flags.
inline void set_gluniform_texture(const glprogram& prog, int var, int varon,
    const gltexture_info& tinfo, uint tunit) {
    set_gluniform_texture(prog, var, tinfo, tunit);
    set_gluniform(prog, varon, is_gltexture_valid(tinfo.txt));
}

// Set uniform texture.
inline void set_gluniform_texture(const glprogram& prog, const std::string& var,
    const gltexture_info& tinfo, uint tunit) {
    auto loc = get_gluniform_location(prog, var);
    if (loc < 0) throw std::runtime_error("bad OpenGL id");
    return set_gluniform_texture(prog, loc, tinfo, tunit);
}
// Set uniform texture with an additionasl texture enable flags.
inline void set_gluniform_texture(const glprogram& prog, const std::string& var,
    const std::string& varon, const gltexture_info& tinfo, uint tunit) {
    auto loc = get_gluniform_location(prog, var);
    if (loc < 0) throw std::runtime_error("bad OpenGL id");
    auto locon = get_gluniform_location(prog, varon);
    if (locon < 0) throw std::runtime_error("bad OpenGL id");
    return set_gluniform_texture(prog, loc, locon, tinfo, tunit);
}

// Binds a buffer to a vertex attribute, or a constant if the buffer is empty.
void set_glattribute(
    const glprogram& prog, int var, const glvertex_buffer& buf, float def);
void set_glattribute(const glprogram& prog, int var, const glvertex_buffer& buf,
    const vec2f& def);
void set_glattribute(const glprogram& prog, int var, const glvertex_buffer& buf,
    const vec3f& def);
void set_glattribute(const glprogram& prog, int var, const glvertex_buffer& buf,
    const vec4f& def);

// Binds a buffer or constant to a vertex attribute.
template <typename T>
inline void set_glattribute(const glprogram& prog, const std::string& var,
    const glvertex_buffer& buf, const T& def) {
    auto loc = get_glattrib_location(prog, var);
    if (loc < 0) throw std::runtime_error("bad OpenGL id");
    set_glattribute(prog, loc, buf, def);
}

// Check whether the program is valid.
inline bool is_glprogram_valid(const glprogram& prog) { return (bool)prog.pid; }

// Binds/unbinds a program.
void bind_glprogram(const glprogram& prog);
void unbind_glprogram(const glprogram& prog);

// Clears OpenGL state.
void clear_program(glprogram& prog);

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMAGE SHADER
// -----------------------------------------------------------------------------
namespace ygl {

// A shader for displaying images.  Members are not part of the public API.
struct glimage_program {
    glprogram prog;       // program
    glvertex_buffer vbo;  // vertex array
    glvertex_buffer ebo;  // element array
};

// Initialize a stdimage program.
glimage_program make_glimage_program();

// Draws an image texture the stdimage program.
void draw_glimage(const glimage_program& prog, const gltexture& txt,
    const vec2i& win_size, const vec2f& offset, float zoom, float exposure = 0,
    float gamma = 1);

// Computes the image uv coordinates corresponding to the view parameters.
inline vec2i get_glimage_coords(
    const vec2f& mouse_pos, const vec2f& offset, float zoom) {
    auto xy = (mouse_pos - offset) / zoom;
    return {(int)round(xy.x), (int)round(xy.y)};
}

// Program to shade surfaces with a physically-based standard shader based on
// Phong/GGX. Members are not part of public API.
struct glsurface_program {
    glprogram prog;  // program
    // uniform variable location
    int eyelight_id, exposure_id, gamma_id, cam_xform_id, cam_xform_inv_id,
        cam_proj_id, lamb_id, lnum_id, lpos_id[16], lke_id[16], ltype_id[16],
        shp_xform_id, shp_normal_offset_id, highlight_id, mtype_id, ke_id,
        kd_id, ks_id, rs_id, op_id, ke_txt_id, ke_txt_on_id, kd_txt_id,
        kd_txt_on_id, ks_txt_id, ks_txt_on_id, rs_txt_id, rs_txt_on_id,
        norm_txt_id, norm_txt_on_id, occ_txt_id, occ_txt_on_id, norm_scale_id,
        occ_scale_id, use_phong_id, double_sided_id, alpha_cutout_id, etype_id,
        efaceted_id;
    // vertex attribute locations
    int pos_id, norm_id, texcoord_id, color_id, tangsp_id;
};

// Initialize a stdsurface shader.
glsurface_program make_glsurface_program();

// Check if the program is valid.
inline bool is_glprogram_valid(const glsurface_program& prog) {
    return is_glprogram_valid(prog.prog);
}

// Starts a frame by setting exposure/gamma values, camera transforms and
// projection. Sets also whether to use full shading or a quick eye light
// preview.
void begin_glsurface_frame(const glsurface_program& prog,
    const mat4f& camera_xform, const mat4f& camera_xform_inv,
    const mat4f& camera_proj, bool shade_eyelight, float exposure = 0,
    float gamma = 2.2f);

// Ends a frame.
void end_glsurface_frame(const glsurface_program& prog);

// Set shading lights and ambient.
void set_glsurface_lights(
    const glsurface_program& prog, const vec3f& amb, const gllights& lights);

// Begins drawing a shape with transform `xform`.
void begin_glsurface_shape(
    const glsurface_program& prog, const mat4f& xform, float normal_offset = 0);

// End shade drawing.
void end_glsurface_shape(const glsurface_program& prog);

// Sets normal offset.
void set_glsurface_normaloffset(
    const glsurface_program& prog, float normal_offset);

// Set the object as highlighted.
void set_glsurface_highlight(
    const glsurface_program& prog, const vec4f& highlight);

// Set material values with emission `ke`, diffuse `kd`, specular `ks` and
// specular roughness `rs`, opacity `op`. Indicates textures ids with the
// correspoinding `XXX_txt` variables. Sets also normal and occlusion
// maps. Works for points/lines/triangles indicated by `etype`, (diffuse for
// points, Kajiya-Kay for lines, GGX/Phong for triangles). Material `type`
// matches the scene material type.
// Mtype: (1) specular-roughness (2) metallic-roughness (3) specular-glossiness
void set_glsurface_material(const glsurface_program& prog, int mtype,
    const vec3f& ke, const vec3f& kd, const vec3f& ks, float rs, float op,
    const gltexture_info& ke_txt, const gltexture_info& kd_txt,
    const gltexture_info& ks_txt, const gltexture_info& rs_txt,
    const gltexture_info& norm_txt, const gltexture_info& occ_txt,
    bool use_phong, bool double_sided, bool alpha_cutout);

// Set constant material with emission `ke` and opacity `op`.
void set_glsurface_constmaterial(
    const glsurface_program& prog, const vec3f& ke, float op);

// Set element properties.
void set_glsurface_elems(
    const glsurface_program& prog, glelem_type etype, bool faceted);

// Set vertex data with buffers for position pos, normals norm, texture
// coordinates texcoord, per-vertex color color and tangent space tangsp.
void set_glsurface_vert(const glsurface_program& prog,
    const glvertex_buffer& pos, const glvertex_buffer& norm,
    const glvertex_buffer& texcoord, const glvertex_buffer& color,
    const glvertex_buffer& tangsp);

}  // namespace ygl

// Forward declaration
struct GLFWwindow;

// -----------------------------------------------------------------------------
// OPENGL WINDOWS
// -----------------------------------------------------------------------------
namespace ygl {

// Forward declaration
struct glwindow;

// Callbacks.
using text_glcallback = std::function<void(unsigned int key)>;
using mouse_glcallback = std::function<void(int button, bool press, int mods)>;
using refresh_glcallback = std::function<void()>;

// OpenGL window. Members are not part of the public API.
struct glwindow {
    GLFWwindow* gwin = nullptr;               // GLFW window
    bool widget_enabled = false;              // whether we have widgets
    text_glcallback text_cb = nullptr;        // text callback
    mouse_glcallback mouse_cb = nullptr;      // mouse callback
    refresh_glcallback refresh_cb = nullptr;  // refresh callback

    ~glwindow();  // cleaup
};

// Initialize a window.
glwindow* make_glwindow(
    int width, int height, const std::string& title, bool opengl4 = true);

// Set window callbacks.
void set_glwindow_callbacks(glwindow* win, text_glcallback text_cb,
    mouse_glcallback mouse_cb, refresh_glcallback refresh_cb);

// Set window title.
void set_glwindow_title(glwindow* win, const std::string& title);

// Event processing.
void wait_glwindow_events(glwindow* win);
void poll_glwindow_events(glwindow* win);
void swap_glwindow_buffers(glwindow* win);
#ifdef __APPLE__
void wait_glwindow_events_timeout(glwindow* win, double timeout_sec);
void post_glwindow_event(glwindow* win);
#endif
// Whether the window should exit the event processing loop.
bool should_glwindow_close(glwindow* win);

// Window/framebuffer size.
vec2i get_glwindow_size(glwindow* win);
vec2i get_glframebuffer_size(glwindow* win);

// Mouse/keyboard state queries.
int get_glmouse_button(glwindow* win);
vec2i get_glmouse_pos(glwindow* win);
vec2f get_glmouse_posf(glwindow* win);
bool get_glkey(glwindow* win, int key);
bool get_glalt_key(glwindow* win);
bool get_glctrl_key(glwindow* win);
bool get_glshift_key(glwindow* win);

// Read pixels.
void take_glscreenshot4b(glwindow* win, int& width, int& height,
    std::vector<vec4b>& img, bool flipy = true, bool back = false);

}  // namespace ygl

// -----------------------------------------------------------------------------
// OPENGL WIDGETS
// -----------------------------------------------------------------------------
namespace ygl {

// Initialize widgets.
void init_imgui(
    glwindow* win, bool light_style = false, bool extra_font = true);

// Begin/end draw widgets.
bool begin_imgui_frame(glwindow* win, const std::string& title);
void end_imgui_frame(glwindow* win);

// Whether widgets are active.
bool get_imgui_active(glwindow* win);

// Horizontal separator.
void draw_imgui_separator(glwindow* win);

// Indent and line continuation widget.
void begin_imgui_indent(glwindow* win);
void end_imgui_indent(glwindow* win);
void continue_imgui_line(glwindow* win);

// Label widget.
void draw_imgui_label(
    glwindow* win, const std::string& lbl, const std::string& msg);

// Checkbox widget
bool draw_imgui_checkbox(glwindow* win, const std::string& lbl, bool& val);
// Text widget.
bool draw_imgui_text(glwindow* win, const std::string& lbl, std::string& str);
bool draw_imgui_multiline_text(
    glwindow* win, const std::string& lbl, std::string& str);

// Drag widget scale (defaults to 1/100).
void draw_drag_speedscale(float scale);
// Drag widget.
bool draw_imgui_dragbox(
    glwindow* win, const std::string& lbl, int& val, int min = 0, int max = 1);
bool draw_imgui_dragbox(glwindow* win, const std::string& lbl, vec2i& val,
    int min = 0, int max = 1);
bool draw_imgui_dragbox(glwindow* win, const std::string& lbl, vec3i& val,
    int min = 0, int max = 1);
bool draw_imgui_dragbox(glwindow* win, const std::string& lbl, vec4i& val,
    int min = 0, int max = 1);
bool draw_imgui_dragbox(glwindow* win, const std::string& lbl, float& val,
    float min = 0, float max = 1);
bool draw_imgui_dragbox(glwindow* win, const std::string& lbl, vec2f& val,
    float min = 0, float max = 1);
bool draw_imgui_dragbox(glwindow* win, const std::string& lbl, vec3f& val,
    float min = 0, float max = 1);
bool draw_imgui_dragbox(glwindow* win, const std::string& lbl, vec4f& val,
    float min = 0, float max = 1);
bool draw_imgui_dragbox(glwindow* win, const std::string& lbl, mat4f& val,
    float min = -1, float max = 1);
bool draw_imgui_dragbox(glwindow* win, const std::string& lbl, frame3f& val,
    float min = -10, float max = 10);

// Color widget.
bool draw_imgui_colorbox(glwindow* win, const std::string& lbl, vec4b& val);
bool draw_imgui_colorbox(glwindow* win, const std::string& lbl, vec4f& val);
bool draw_imgui_colorbox(glwindow* win, const std::string& lbl, vec3f& val);
bool draw_hdr_color_widget(
    glwindow* win, const std::string& lbl, vec3f& val, float max = 10);

// Combo widget.
bool begin_imgui_combobox(
    glwindow* win, const std::string& lbl, const std::string& label);
bool draw_imgui_item(
    glwindow* win, const std::string& label, int idx, bool selected);
void end_imgui_combobox(glwindow* win);
template <typename T>
bool draw_imgui_item(
    glwindow* win, const std::string& label, int idx, T& val, const T& item) {
    auto selected = draw_imgui_item(win, label, idx, val == item);
    if (selected) val = item;
    return selected;
}
inline bool draw_imgui_combobox(glwindow* win, const std::string& lbl,
    std::string& val, const std::vector<std::string>& labels);
template <typename T>
inline bool draw_imgui_combobox(glwindow* win, const std::string& lbl, T& val,
    const std::map<T, std::string>& labels);
template <typename T>
inline bool draw_imgui_combobox(glwindow* win, const std::string& lbl, T*& val,
    const std::vector<T*>& vals, bool extra = true, T* extra_val = nullptr);

// Button widget.
bool draw_imgui_button(glwindow* win, const std::string& lbl);

// Collapsible header widget.
bool draw_imgui_header(glwindow* win, const std::string& lbl);

// Tree widget.
bool begin_imgui_tree(glwindow* win, const std::string& lbl);
void end_imgui_tree(glwindow* win);
bool begin_imgui_tree(
    glwindow* win, const std::string& lbl, void*& selection, void* content);
template <typename T>
inline bool begin_imgui_tree(
    glwindow* win, const std::string& lbl, T*& selection, T* content) {
    auto sel = selection;
    auto open = begin_imgui_tree(win, lbl, (void*&)sel, (void*)content);
    if (sel == content) selection = content;
    return open;
}
void end_imgui_tree(glwindow* win, void* content);
void draw_imgui_tree_leaf(
    glwindow* win, const std::string& lbl, void*& selection, void* content);
template <typename T>
inline void draw_imgui_tree_leaf(
    glwindow* win, const std::string& lbl, void*& selection, T* content) {
    auto sel = selection;
    draw_imgui_tree_leaf(win, lbl, sel, content);
    if (sel == content) selection = content;
}
void draw_imgui_tree_leaf(glwindow* win, const std::string& lbl,
    void*& selection, void* content, const vec4f& col);

// Image widget.
void draw_imgui_imagebox(
    glwindow* win, int tid, const vec2i& size, const vec2i& imsize);
void draw_imgui_imagebox(glwindow* win, gltexture& txt, const vec2i& size);

// Scroll region widget.
void begin_imgui_scrollarea(
    glwindow* win, const std::string& lbl, int height, bool border);
void end_imgui_scrollarea(glwindow* win);
void move_imgui_scrollarea(glwindow* win);

// Group ids widget.
void push_imgui_groupid(glwindow* win, int gid);
void push_imgui_groupid(glwindow* win, void* gid);
void push_imgui_groupid(glwindow* win, const void* gid);
void push_imgui_groupid(glwindow* win, const char* gid);
void pop_imgui_groupid(glwindow* win);

// Widget style.
void push_imgui_style(glwindow* win, const vec4f& color);
void pop_imgui_style(glwindow* win);

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR OPENGL WIDGETS
// -----------------------------------------------------------------------------
namespace ygl {

// Combo widget.
inline bool draw_imgui_combobox(glwindow* win, const std::string& lbl,
    std::string& val, const std::vector<std::string>& labels) {
    if (!begin_imgui_combobox(win, lbl, val)) return false;
    auto old_val = val;
    for (auto i = 0; i < labels.size(); i++) {
        draw_imgui_item(win, labels[i], i, val, labels[i]);
    }
    end_imgui_combobox(win);
    return val != old_val;
}

// Combo widget.
template <typename T>
inline bool draw_imgui_combobox(glwindow* win, const std::string& lbl, T& val,
    const std::map<T, std::string>& labels) {
    if (!begin_imgui_combobox(win, lbl, labels.at(val))) return false;
    auto old_val = val;
    auto lid = 0;
    for (auto& kv : labels)
        draw_imgui_item(win, kv.second, lid++, val, kv.first);
    end_imgui_combobox(win);
    return val != old_val;
}

// Combo widget
template <typename T>
inline bool draw_imgui_combobox(glwindow* win, const std::string& lbl, T*& val,
    const std::vector<T*>& vals, bool extra, T* extra_val) {
    if (!begin_imgui_combobox(win, lbl, (val) ? val->name : "<none>"))
        return false;
    auto old_val = val;
    if (extra)
        draw_imgui_item(
            win, (extra_val) ? extra_val->name : "<none>", -1, val, extra_val);
    for (auto i = 0; i < vals.size(); i++) {
        draw_imgui_item(win, vals[i]->name, i, val, vals[i]);
    }
    end_imgui_combobox(win);
    return val != old_val;
}

}  // namespace ygl

#endif
