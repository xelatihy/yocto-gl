//
// Yocto/OpenGL: Utilities to use OpenGL 3, GLFW and ImGui.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2019 Fabio Pellacini
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

#ifndef _YOCTO_OPENGL_
#define _YOCTO_OPENGL_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <functional>
#include <memory>

#include "../yocto/yocto_image.h"
#include "../yocto/yocto_math.h"

// forward declaration
struct GLFWwindow;

// -----------------------------------------------------------------------------
// IMAGE DRAWING
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::unique_ptr;

// OpenGL image data
struct opengl_image {
  vec2i size() const { return texture_size; }
        operator bool() const { return (bool)texture_id; }

  uint  program_id     = 0;
  uint  vertex_id      = 0;
  uint  fragment_id    = 0;
  uint  array_id       = 0;
  uint  texcoords_id   = 0;
  uint  triangles_id   = 0;
  uint  texture_id     = 0;
  vec2i texture_size   = {0, 0};
  bool  texture_linear = false;
  bool  texture_mipmap = false;

  ~opengl_image();
};

// create image drawing program
opengl_image* make_glimage();

// update image data
void set_glimage(opengl_image* glimage, const image<vec4f>& img,
    bool linear = false, bool mipmap = false);
void set_glimage(opengl_image* glimage, const image<vec4b>& img,
    bool linear = false, bool mipmap = false);

// OpenGL image drawing params
struct draw_glimage_params {
  vec2i window      = {512, 512};
  vec4i framebuffer = {0, 0, 512, 512};
  vec2f center      = {0, 0};
  float scale       = 1;
  bool  fit         = true;
  bool  checker     = true;
  float border_size = 2;
  vec4f background  = {0.15f, 0.15f, 0.15f, 1.0f};
};

// draw image
void draw_glimage(opengl_image* glimage, const draw_glimage_params& params);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE DRAWING
// -----------------------------------------------------------------------------
namespace yocto {

// Opengl caemra
struct opengl_camera {
  frame3f frame  = identity3x4f;
  float   lens   = 0.050;
  float   asepct = 1;
  float   film   = 0.036;
  float   near   = 0.001;
  float   far    = 10000;
};

// Opengl shape
struct opengl_shape {
  int  positions_num = 0;
  uint positions_id  = 0;
  int  normals_num   = 0;
  uint normals_id    = 0;
  int  texcoords_num = 0;
  uint texcoords_id  = 0;
  int  colors_num    = 0;
  uint colors_id     = 0;
  int  tangents_num  = 0;
  uint tangents_id   = 0;
  int  points_num    = 0;
  uint points_id     = 0;
  int  lines_num     = 0;
  uint lines_id      = 0;
  int  triangles_num = 0;
  uint triangles_id  = 0;
  int  quads_num     = 0;
  uint quads_id      = 0;
  int  edges_num     = 0;
  uint edges_id      = 0;
};

// OpenGL texture
struct opengl_texture_ {
  uint  texture_id = 0;
  vec2i size       = {0, 0};
  bool  mipmap     = false;
  bool  linear     = false;
  bool  is_srgb    = false;
  bool  is_float   = false;
};

// Opengl material
struct opengl_material {
  vec3f emission      = zero3f;
  vec3f diffuse       = zero3f;
  vec3f specular      = zero3f;
  float metallic      = 0;
  float roughness     = 0;
  float opacity       = 1;
  int   emission_map  = -1;
  int   diffuse_map   = -1;
  int   specular_map  = -1;
  int   metallic_map  = -1;
  int   roughness_map = -1;
  int   normal_map    = -1;
  bool  gltf_textures = false;
};

// Opengl instance group
struct opengl_instance {
  frame3f frame       = identity3x4f;
  int     shape       = 0;
  int     material    = 0;
  bool    highlighted = false;
};

// Opengl light
struct opengl_light {
  vec3f position = zero3f;
  vec3f emission = zero3f;
  int   type     = 0;
};

// Opengl scene
struct opengl_scene {
  vector<opengl_camera>   _cameras   = {};
  vector<opengl_instance> _instances = {};
  vector<opengl_shape>    _shapes    = {};
  vector<opengl_material> _materials = {};
  vector<opengl_texture_> _textures  = {};
  vector<opengl_light>    _lights    = {};

  // OpenGL state
  uint program_id  = 0;
  uint vertex_id   = 0;
  uint fragment_id = 0;
  uint array_id    = 0;

  // cleanup
  ~opengl_scene();
};

// Draw options
struct draw_glscene_params {
  int   camera           = 0;
  int   resolution       = 1280;
  bool  wireframe        = false;
  bool  edges            = false;
  float edge_offset      = 0.01f;
  bool  eyelight         = false;
  float exposure         = 0;
  float gamma            = 2.2f;
  vec3f ambient          = {0, 0, 0};
  bool  double_sided     = true;
  bool  non_rigid_frames = true;
  float near             = 0.01f;
  float far              = 10000.0f;
  vec4f background       = vec4f{0.15f, 0.15f, 0.15f, 1.0f};
};

// Initialize an OpenGL scene
void make_glscene(opengl_scene& scene);

// add camera
int  add_glcamera(opengl_scene& scene);
void set_glcamera_frame(opengl_scene& scene, int idx, const frame3f frame);
void set_glcamera_lens(
    opengl_scene& scene, int idx, float lens, float asepct, float film);
void set_glcamera_planes(opengl_scene& scene, int idx, float near, float far);
void clear_glcameras(opengl_scene& scene);

// add texture
int  add_gltexture(opengl_scene& scene);
void set_gltexture(
    opengl_scene& scene, int idx, const image<vec4b>& img, bool as_srgb = true);
void set_gltexture(opengl_scene& scene, int idx, const image<vec4f>& img,
    bool as_float = false);
void clear_gltextures(opengl_scene& scene);

// add material
int  add_glmaterial(opengl_scene& scene);
void set_glmaterial_emission(
    opengl_scene& scene, int idx, const vec3f& emission, int emission_txt = -1);
void set_glmaterial_diffuse(
    opengl_scene& scene, int idx, const vec3f& diffuse, int diffuse_txt = -1);
void set_glmaterial_specular(
    opengl_scene& scene, int idx, const vec3f& specular, int specular_txt = -1);
void set_glmaterial_roughness(
    opengl_scene& scene, int idx, float roughness, int roughness_txt = -1);
void set_glmaterial_opacity(
    opengl_scene& scene, int idx, float opacity, int opacity_txt = -1);
void set_glmaterial_metallic(
    opengl_scene& scene, int idx, float metallic, int metallic_txt = -1);
void set_glmaterial_normalmap(opengl_scene& scene, int idx, int normal_txt);
void set_glmaterial_gltftextures(
    opengl_scene& scene, int idx, bool gltf_textures);
void clean_glmaterias(opengl_scene& scene);

// add shape
int  add_glshape(opengl_scene& scene);
void set_glshape_positions(
    opengl_scene& scene, int idx, const vector<vec3f>& positions);
void set_glshape_normals(
    opengl_scene& scene, int idx, const vector<vec3f>& normals);
void set_glshape_texcoords(
    opengl_scene& scene, int idx, const vector<vec2f>& texcoords);
void set_glshape_colors(
    opengl_scene& scene, int idx, const vector<vec4f>& colors);
void set_glshape_tangents(
    opengl_scene& scene, int idx, const vector<vec4f>& tangents);
void set_glshape_points(
    opengl_scene& scene, int idx, const vector<int>& points);
void set_glshape_lines(
    opengl_scene& scene, int idx, const vector<vec2i>& lines);
void set_glshape_triangles(
    opengl_scene& scene, int idx, const vector<vec3i>& triangles);
void set_glshape_quads(
    opengl_scene& scene, int idx, const vector<vec4i>& quads);
void set_glshape_edges(
    opengl_scene& scene, int idx, const vector<vec2i>& edges);
void clean_glshapes(opengl_scene& scene);

// add instance
int  add_glinstance(opengl_scene& scene);
void set_glinstance_frame(opengl_scene& scene, int idx, const frame3f& frame);
void set_glinstance_shape(opengl_scene& scene, int idx, int shape);
void set_glinstance_material(opengl_scene& scene, int idx, int material);
void clear_glinstances(opengl_scene& scene);

// add light
int  add_gllight(opengl_scene& scene, const vec3f& position,
     const vec3f& emission, bool directional);
void set_gllight(opengl_scene& scene, int idx, const vec3f& position,
    const vec3f& emission, bool directional);
void clear_gllights(opengl_scene& scene);
bool has_max_gllights(opengl_scene& scene);

// Draw an OpenGL scene
void draw_glscene(opengl_scene& state, const vec4i& viewport,
    const draw_glscene_params& params);

}  // namespace yocto

// -----------------------------------------------------------------------------
// OPENGL WINDOW
// -----------------------------------------------------------------------------
namespace yocto {

// Forward declaration of OpenGL window
struct opengl_window;

// Input state
struct opengl_input {
  bool     mouse_left     = false;  // left button
  bool     mouse_right    = false;  // right button
  bool     mouse_middle   = false;  // middle button
  vec2f    mouse_pos      = {};     // position excluding widgets
  vec2f    mouse_last     = {};     // last mouse position excluding widgets
  vec2f    mouse_delta    = {};     // last mouse delta excluding widgets
  bool     modifier_alt   = false;  // alt modifier
  bool     modifier_ctrl  = false;  // ctrl modifier
  bool     modifier_shift = false;  // shift modifier
  bool     widgets_active = false;  // widgets are active
  uint64_t clock_now      = 0;      // clock now
  uint64_t clock_last     = 0;      // clock last
  double   time_now       = 0;      // time now
  double   time_delta     = 0;      // time delta
};

// Draw callback called every frame and when resizing
using draw_glcallback =
    std::function<void(const opengl_window&, vec2i window, vec4i viewport)>;
// Draw callback for drawing widgets
using widgets_glcallback = std::function<void(const opengl_window&)>;
// Drop callback that returns that list of dropped strings.
using drop_glcallback =
    std::function<void(const opengl_window&, const vector<string>&)>;
// Key callback that returns ASCII key, pressed/released flag and modifier keys
using key_glcallback =
    std::function<void(const opengl_window& win, int key, bool pressed)>;
// Mouse click callback that returns left/right button, pressed/released flag,
// modifier keys
using click_glcallback =
    std::function<void(const opengl_window&, bool left, bool pressed)>;
// Scroll callback that returns scroll amount
using scroll_glcallback =
    std::function<void(const opengl_window&, float amount)>;
// Update functions called every frame
using uiupdate_glcallback =
    std::function<void(const opengl_window&, const opengl_input& input)>;
// Update functions called every frame
using update_glcallback = std::function<void(const opengl_window&)>;

// OpenGL window wrapper
struct opengl_window {
  GLFWwindow*         win           = nullptr;
  draw_glcallback     draw_cb       = {};
  widgets_glcallback  widgets_cb    = {};
  drop_glcallback     drop_cb       = {};
  key_glcallback      key_cb        = {};
  click_glcallback    click_cb      = {};
  scroll_glcallback   scroll_cb     = {};
  update_glcallback   update_cb     = {};
  uiupdate_glcallback uiupdate_cb   = {};
  int                 widgets_width = 0;
  bool                widgets_left  = true;
  opengl_input        input         = {};
  vec4f               background    = {0.15f, 0.15f, 0.15f, 1.0f};
};

// Windows initialization
void init_glwindow(opengl_window& win, const vec2i& size, const string& title);

// Window cleanup
void delete_glwindow(opengl_window& win);

// Set callbacks
void set_draw_glcallback(opengl_window& win, draw_glcallback draw_cb);
void set_widgets_glcallback(opengl_window& win, widgets_glcallback widgets_cb);
void set_drop_glcallback(opengl_window& win, drop_glcallback drop_cb);
void set_key_glcallback(opengl_window& win, key_glcallback cb);
void set_click_glcallback(opengl_window& win, click_glcallback cb);
void set_scroll_glcallback(opengl_window& win, scroll_glcallback cb);
void set_uiupdate_glcallback(opengl_window& win, uiupdate_glcallback cb);
void set_update_glcallback(opengl_window& win, update_glcallback cb);

// Run loop
void run_ui(opengl_window& win);

vec2i get_glwindow_size(const opengl_window& win, bool ignore_widgets = true);
vec2i get_glframebuffer_size(
    const opengl_window& win, bool ignore_widgets = true);
vec4i get_glframebuffer_viewport(
    const opengl_window& win, bool ignore_widgets = true);

bool should_glwindow_close(const opengl_window& win);
void set_glwindow_close(const opengl_window& win, bool close);

vec2f get_glmouse_pos(const opengl_window& win, bool ignore_widgets = true);
vec2f get_glmouse_pos_normalized(
    const opengl_window& win, bool ignore_widgets = true);

bool get_glmouse_left(const opengl_window& win);
bool get_glmouse_right(const opengl_window& win);
bool get_glalt_key(const opengl_window& win);
bool get_glshift_key(const opengl_window& win);
bool get_glctrl_key(const opengl_window& win);

void process_glevents(const opengl_window& win, bool wait = false);
void swap_glbuffers(const opengl_window& win);

}  // namespace yocto

// -----------------------------------------------------------------------------
// OPENGL WIDGETS
// -----------------------------------------------------------------------------
namespace yocto {

void init_glwidgets(opengl_window& win, int width = 320, bool left = true);
bool get_glwidgets_active(const opengl_window& win);

void begin_glwidgets(const opengl_window& win);
void end_glwidgets(const opengl_window& win);

bool begin_glwidgets_window(const opengl_window& win, const char* title);

bool begin_glheader(const opengl_window& win, const char* title);
void end_glheader(const opengl_window& win);

void open_glmodal(const opengl_window& win, const char* lbl);
void clear_glmodal(const opengl_window& win);
bool begin_glmodal(const opengl_window& win, const char* lbl);
void end_glmodal(const opengl_window& win);
bool is_glmodal_open(const opengl_window& win, const char* lbl);

bool draw_glmessages(const opengl_window& win);
void push_glmessage(const string& message);
void push_glmessage(const opengl_window& win, const string& message);
bool draw_glfiledialog(const opengl_window& win, const char* lbl, string& path,
    bool save, const string& dirname, const string& filename,
    const string& filter);
bool draw_glfiledialog_button(const opengl_window& win, const char* button_lbl,
    bool button_active, const char* lbl, string& path, bool save,
    const string& dirname, const string& filename, const string& filter);

void draw_gllabel(
    const opengl_window& win, const char* lbl, const string& text);

bool begin_header_widget(const opengl_window& win, const char* label);
void end_header_widget(const opengl_window& win);

void draw_glseparator(const opengl_window& win);
void continue_glline(const opengl_window& win);

bool draw_glbutton(const opengl_window& win, const char* lbl);
bool draw_glbutton(const opengl_window& win, const char* lbl, bool enabled);

bool draw_gltextinput(const opengl_window& win, const char* lbl, string& value);
bool draw_glslider(const opengl_window& win, const char* lbl, float& value,
    float min, float max);
bool draw_glslider(const opengl_window& win, const char* lbl, vec2f& value,
    float min, float max);
bool draw_glslider(const opengl_window& win, const char* lbl, vec3f& value,
    float min, float max);
bool draw_glslider(const opengl_window& win, const char* lbl, vec4f& value,
    float min, float max);

bool draw_glslider(
    const opengl_window& win, const char* lbl, int& value, int min, int max);
bool draw_glslider(
    const opengl_window& win, const char* lbl, vec2i& value, int min, int max);
bool draw_glslider(
    const opengl_window& win, const char* lbl, vec3i& value, int min, int max);
bool draw_glslider(
    const opengl_window& win, const char* lbl, vec4i& value, int min, int max);

bool draw_gldragger(const opengl_window& win, const char* lbl, float& value,
    float speed = 1.0f, float min = 0.0f, float max = 0.0f);
bool draw_gldragger(const opengl_window& win, const char* lbl, vec2f& value,
    float speed = 1.0f, float min = 0.0f, float max = 0.0f);
bool draw_gldragger(const opengl_window& win, const char* lbl, vec3f& value,
    float speed = 1.0f, float min = 0.0f, float max = 0.0f);
bool draw_gldragger(const opengl_window& win, const char* lbl, vec4f& value,
    float speed = 1.0f, float min = 0.0f, float max = 0.0f);

bool draw_gldragger(const opengl_window& win, const char* lbl, int& value,
    float speed = 1, int min = 0, int max = 0);
bool draw_gldragger(const opengl_window& win, const char* lbl, vec2i& value,
    float speed = 1, int min = 0, int max = 0);
bool draw_gldragger(const opengl_window& win, const char* lbl, vec3i& value,
    float speed = 1, int min = 0, int max = 0);
bool draw_gldragger(const opengl_window& win, const char* lbl, vec4i& value,
    float speed = 1, int min = 0, int max = 0);

bool draw_glcheckbox(const opengl_window& win, const char* lbl, bool& value);

bool draw_glcoloredit(const opengl_window& win, const char* lbl, vec3f& value);
bool draw_glcoloredit(const opengl_window& win, const char* lbl, vec4f& value);

bool draw_glhdrcoloredit(
    const opengl_window& win, const char* lbl, vec3f& value);
bool draw_glhdrcoloredit(
    const opengl_window& win, const char* lbl, vec4f& value);

bool draw_glcombobox(const opengl_window& win, const char* lbl, int& idx,
    const vector<string>& labels);
bool draw_glcombobox(const opengl_window& win, const char* lbl, string& value,
    const vector<string>& labels);
bool draw_glcombobox(const opengl_window& win, const char* lbl, int& idx,
    int num, const std::function<const char*(int)>& labels,
    bool include_null = false);

template <typename T>
inline bool draw_glcombobox(const opengl_window& win, const char* lbl, int& idx,
    const vector<T>& vals, bool include_null = false) {
  return draw_glcombobox(win, lbl, idx, (int)vals.size(),
      [&](int idx) { return vals[idx].name.c_str(); }, include_null);
}

void draw_glhistogram(
    const opengl_window& win, const char* lbl, const vector<float>& values);
void draw_glhistogram(
    const opengl_window& win, const char* lbl, const vector<vec2f>& values);
void draw_glhistogram(
    const opengl_window& win, const char* lbl, const vector<vec3f>& values);
void draw_glhistogram(
    const opengl_window& win, const char* lbl, const vector<vec4f>& values);

void log_glinfo(const opengl_window& win, const string& msg);
void log_glerror(const opengl_window& win, const string& msg);
void clear_gllogs(const opengl_window& win);
void draw_gllog(const opengl_window& win);

}  // namespace yocto

#endif
