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
  opengl_image() {}
  opengl_image(const opengl_image&) = delete;
  opengl_image& operator=(const opengl_image&) = delete;

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
void init_glimage(opengl_image& glimage);
bool is_initialized(const opengl_image& glimage);

// update image data
void set_glimage(opengl_image& glimage, const image<vec4f>& img,
    bool linear = false, bool mipmap = false);
void set_glimage(opengl_image& glimage, const image<vec4b>& img,
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
void draw_glimage(opengl_image& glimage, const draw_glimage_params& params);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE DRAWING
// -----------------------------------------------------------------------------
namespace yocto {

// Opengl caemra
struct opengl_camera {
  frame3f frame  = identity3x4f;
  float   lens   = 0.050;
  float   aspect = 1.000;
  float   film   = 0.036;
  float   near   = 0.001;
  float   far    = 10000;
};

// OpenGL texture
struct opengl_texture {
  uint  texture_id = 0;
  vec2i size       = {0, 0};
  bool  is_srgb    = false;
  bool  is_float   = false;

  opengl_texture() {}
  opengl_texture(const opengl_texture&) = delete;
  opengl_texture& operator=(opengl_texture&) = delete;
  ~opengl_texture();
  opengl_texture(opengl_texture&&);
  opengl_texture& operator=(opengl_texture&&);
};

// Opengl material
struct opengl_material {
  vec3f emission      = {0, 0, 0};
  vec3f diffuse       = {0, 0, 0};
  vec3f specular      = {0, 0, 0};
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

// Opengl shape
struct opengl_shape {
  // shape properties
  frame3f         frame       = identity3x4f;
  vector<frame3f> instances   = {};
  opengl_material material    = {};
  bool            hidden      = false;
  bool            highlighted = false;

  // vertex buffers
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

  opengl_shape() {}
  opengl_shape(const opengl_shape&) = delete;
  opengl_shape& operator=(const opengl_shape&) = delete;
  ~opengl_shape();
  opengl_shape(opengl_shape&&);
  opengl_shape& operator=(opengl_shape&&);
};

// Opengl light
struct opengl_light {
  vec3f position = {0, 0, 0};
  vec3f emission = {0, 0, 0};
  int   type     = 0;
};

// Opengl scene
struct opengl_scene {
  opengl_scene() {}
  opengl_scene(const opengl_scene&) = delete;
  opengl_scene& operator=(const opengl_scene&) = delete;
  ~opengl_scene();

  vector<opengl_camera>  _cameras  = {};
  vector<opengl_shape>   _shapes   = {};
  vector<opengl_texture> _textures = {};
  vector<opengl_light>   _lights   = {};

  // OpenGL state
  uint program_id  = 0;
  uint vertex_id   = 0;
  uint fragment_id = 0;
  uint array_id    = 0;
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
void init_glscene(opengl_scene& scene);
bool is_initialized(opengl_scene& scene);

// add camera
int  add_camera(opengl_scene& scene);
void set_camera_frame(opengl_scene& scene, int idx, const frame3f& frame);
void set_camera_lens(opengl_scene& scene, int idx, float lens);
void set_camera_aspect(opengl_scene& scene, int idx, float aspect);
void set_camera_film(opengl_scene& scene, int idx, float film);
void set_camera_near(opengl_scene& scene, int idx, float near);
void set_camera_far(opengl_scene& scene, int idx, float far);
void clear_cameras(opengl_scene& scene);

// add texture
int  add_texture(opengl_scene& scene);
void set_texture(
    opengl_scene& scene, int idx, const image<vec4b>& img, bool as_srgb = true);
void set_texture(opengl_scene& scene, int idx, const image<vec4f>& img,
    bool as_float = false);
void clear_textures(opengl_scene& scene);

// add shape
int  add_shape(opengl_scene& scene);
void set_shape_points(opengl_scene& scene, int idx, const vector<int>& points);
void set_shape_lines(opengl_scene& scene, int idx, const vector<vec2i>& lines);
void set_shape_triangles(
    opengl_scene& scene, int idx, const vector<vec3i>& triangles);
void set_shape_quads(opengl_scene& scene, int idx, const vector<vec4i>& quads);
void set_shape_positions(
    opengl_scene& scene, int idx, const vector<vec3f>& positions);
void set_shape_normals(
    opengl_scene& scene, int idx, const vector<vec3f>& normals);
void set_shape_texcoords(
    opengl_scene& scene, int idx, const vector<vec2f>& texcoords);
void set_shape_colors(
    opengl_scene& scene, int idx, const vector<vec4f>& colors);
void set_shape_tangents(
    opengl_scene& scene, int idx, const vector<vec4f>& tangents);
void set_shape_frame(opengl_scene& scene, int idx, const frame3f& frame);
void set_shape_instances(
    opengl_scene& scene, int idx, const vector<frame3f>& instances);
void set_material_emission(
    opengl_scene& scene, int idx, const vec3f& emission, int emission_txt = -1);
void set_material_diffuse(
    opengl_scene& scene, int idx, const vec3f& diffuse, int diffuse_txt = -1);
void set_material_specular(
    opengl_scene& scene, int idx, const vec3f& specular, int specular_txt = -1);
void set_material_roughness(
    opengl_scene& scene, int idx, float roughness, int roughness_txt = -1);
void set_material_opacity(
    opengl_scene& scene, int idx, float opacity, int opacity_txt = -1);
void set_material_metallic(
    opengl_scene& scene, int idx, float metallic, int metallic_txt = -1);
void set_material_normalmap(opengl_scene& scene, int idx, int normal_txt);
void set_material_gltftextures(
    opengl_scene& scene, int idx, bool gltf_textures);
void set_shape_hidden(opengl_scene& scene, int idx, bool hidden);
void set_shape_highlighted(opengl_scene& scene, int idx, bool highlighted);
void clear_shapes(opengl_scene& scene);

// add light
int  add_light(opengl_scene& scene);
void set_light_position(opengl_scene& scene, int idx, const vec3f& position);
void set_light_emission(opengl_scene& scene, int idx, const vec3f& emission);
void set_light_directional(opengl_scene& scene, int idx, bool directional);
void set_light(opengl_scene& scene, int idx, const vec3f& position,
    const vec3f& emission, bool directional);
void clear_lights(opengl_scene& scene);
bool has_max_lights(opengl_scene& scene);

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
  bool     mouse_left           = false;  // left button
  bool     mouse_right          = false;  // right button
  bool     mouse_middle         = false;  // middle button
  vec2f    mouse_pos            = {};     // position excluding widgets
  vec2f    mouse_last           = {};  // last mouse position excluding widgets
  vec2f    mouse_delta          = {};  // last mouse delta excluding widgets
  bool     modifier_alt         = false;         // alt modifier
  bool     modifier_ctrl        = false;         // ctrl modifier
  bool     modifier_shift       = false;         // shift modifier
  bool     widgets_active       = false;         // widgets are active
  uint64_t clock_now            = 0;             // clock now
  uint64_t clock_last           = 0;             // clock last
  double   time_now             = 0;             // time now
  double   time_delta           = 0;             // time delta
  vec2i    window_size          = {0, 0};        // window size
  vec4i    framebuffer_viewport = {0, 0, 0, 0};  // framebuffer viewport
};

// Draw callback called every frame and when resizing
using draw_glcallback =
    std::function<void(const opengl_window&, const opengl_input& input)>;
// Draw callback for drawing widgets
using widgets_glcallback =
    std::function<void(const opengl_window&, const opengl_input& input)>;
// Drop callback that returns that list of dropped strings.
using drop_glcallback = std::function<void(
    const opengl_window&, const vector<string>&, const opengl_input& input)>;
// Key callback that returns ASCII key, pressed/released flag and modifier keys
using key_glcallback = std::function<void(
    const opengl_window&, int key, bool pressed, const opengl_input& input)>;
// Mouse click callback that returns left/right button, pressed/released flag,
// modifier keys
using click_glcallback = std::function<void(
    const opengl_window&, bool left, bool pressed, const opengl_input& input)>;
// Scroll callback that returns scroll amount
using scroll_glcallback = std::function<void(
    const opengl_window&, float amount, const opengl_input& input)>;
// Update functions called every frame
using uiupdate_glcallback =
    std::function<void(const opengl_window&, const opengl_input& input)>;
// Update functions called every frame
using update_glcallback =
    std::function<void(const opengl_window&, const opengl_input& input)>;

// OpenGL window wrapper
struct opengl_window {
  GLFWwindow*         win           = nullptr;
  string              title         = "";
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
void init_glwindow(opengl_window& win, const vec2i& size, const string& title,
    bool widgets, int widgets_width = 320, bool widgets_left = true);

// Window cleanup
void clear_glwindow(opengl_window& win);

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
void set_close(const opengl_window& win, bool close);

}  // namespace yocto

// -----------------------------------------------------------------------------
// OPENGL WIDGETS
// -----------------------------------------------------------------------------
namespace yocto {

bool begin_glheader(const opengl_window& win, const char* title);
void end_glheader(const opengl_window& win);

void draw_gllabel(
    const opengl_window& win, const char* lbl, const string& text);

void draw_glseparator(const opengl_window& win);
void continue_glline(const opengl_window& win);

bool draw_glbutton(
    const opengl_window& win, const char* lbl, bool enabled = true);

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
  return draw_glcombobox(
      win, lbl, idx, (int)vals.size(),
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

bool draw_glmessages(const opengl_window& win);
void push_glmessage(const opengl_window& win, const string& message);
bool draw_glfiledialog(const opengl_window& win, const char* lbl, string& path,
    bool save, const string& dirname, const string& filename,
    const string& filter);
bool draw_glfiledialog_button(const opengl_window& win, const char* button_lbl,
    bool button_active, const char* lbl, string& path, bool save,
    const string& dirname, const string& filename, const string& filter);

void log_glinfo(const opengl_window& win, const string& msg);
void log_glerror(const opengl_window& win, const string& msg);
void clear_gllogs(const opengl_window& win);
void draw_gllog(const opengl_window& win);

}  // namespace yocto

#endif
