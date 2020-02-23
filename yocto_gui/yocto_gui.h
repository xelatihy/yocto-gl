//
// Yocto/Gui: Utilities for writing simple native GUIs.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2020 Fabio Pellacini
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

#ifndef _YOCTO_GUI_
#define _YOCTO_GUI_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <functional>
#include <memory>
#include <vector>

#include "../yocto_gl/yocto_image.h"
#include "../yocto_gl/yocto_math.h"

// forward declaration
struct GLFWwindow;

// -----------------------------------------------------------------------------
// IMAGE DRAWING
// -----------------------------------------------------------------------------
namespace ygui {

// Math defitions
using ym::bbox3f;
using ym::byte;
using ym::frame3f;
using ym::identity3x4f;
using ym::uint;
using ym::vec2f;
using ym::vec2i;
using ym::vec3b;
using ym::vec3f;
using ym::vec3i;
using ym::vec4b;
using ym::vec4f;
using ym::vec4i;

// OpenGL image data
struct image {
  image() {}
  image(const image&) = delete;
  image& operator=(const image&) = delete;

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

  ~image();
};

// create image drawing program
void init_glimage(ygui::image* image);
bool is_initialized(const ygui::image* image);

// update image data
void set_glimage(ygui::image* image, const yimg::image<vec4f>& img,
    bool linear = false, bool mipmap = false);
void set_glimage(ygui::image* image, const yimg::image<vec4b>& img,
    bool linear = false, bool mipmap = false);

// OpenGL image drawing params
struct image_params {
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
void draw_glimage(ygui::image* image, const image_params& params);

}  // namespace ygui

// -----------------------------------------------------------------------------
// SCENE DRAWING
// -----------------------------------------------------------------------------
namespace ygui {

// Opengl caemra
struct camera {
  frame3f frame  = identity3x4f;
  float   lens   = 0.050;
  float   aspect = 1.000;
  float   film   = 0.036;
  float   near   = 0.001;
  float   far    = 10000;
};

// OpenGL texture
struct texture {
  uint  texture_id = 0;
  vec2i size       = {0, 0};
  int   nchan      = 0;
  bool  is_srgb    = false;
  bool  is_float   = false;

  texture() {}
  texture(const texture&) = delete;
  texture& operator=(texture&) = delete;
  ~texture();
};

// Opengl material
struct material {
  // material
  vec3f          emission      = {0, 0, 0};
  vec3f          color         = {0, 0, 0};
  float          metallic      = 0;
  float          roughness     = 0;
  float          specular      = 0;
  float          opacity       = 1;
  ygui::texture* emission_tex  = nullptr;
  ygui::texture* color_tex     = nullptr;
  ygui::texture* metallic_tex  = nullptr;
  ygui::texture* roughness_tex = nullptr;
  ygui::texture* specular_tex  = nullptr;
  ygui::texture* opacity_tex   = nullptr;
  ygui::texture* normal_tex    = nullptr;
};

// Opengl shape
struct shape {
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

  shape() {}
  shape(const shape&) = delete;
  shape& operator=(const shape&) = delete;
  ~shape();
};

// Opengl instance
struct instance {
  // instancing not supported yet
};

// Opengl object
struct object {
  // object properties
  frame3f         frame       = identity3x4f;
  ygui::shape*    shape       = nullptr;
  ygui::material* material    = nullptr;
  ygui::instance* instance    = nullptr;
  bool            hidden      = false;
  bool            highlighted = false;
};

// Opengl light
struct light {
  vec3f position = {0, 0, 0};
  vec3f emission = {0, 0, 0};
  int   type     = 0;
};

// Opengl scene
struct scene {
  scene() {}
  scene(const scene&) = delete;
  scene& operator=(const scene&) = delete;
  ~scene();

  std::vector<ygui::camera*>   cameras   = {};
  std::vector<ygui::object*>   objects   = {};
  std::vector<ygui::shape*>    shapes    = {};
  std::vector<ygui::material*> materials = {};
  std::vector<ygui::instance*> instances = {};
  std::vector<ygui::texture*>  textures  = {};
  std::vector<ygui::light*>    lights    = {};

  // OpenGL state
  uint program_id  = 0;
  uint vertex_id   = 0;
  uint fragment_id = 0;
  uint array_id    = 0;
};

// Draw options
struct scene_params {
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
void init_glscene(ygui::scene* glscene);
bool is_initialized(const ygui::scene* glscene);

// add scene elements
ygui::camera*   add_camera(ygui::scene* scene);
ygui::texture*  add_texture(ygui::scene* scene);
ygui::material* add_material(ygui::scene* scene);
ygui::shape*    add_shape(ygui::scene* scene);
ygui::instance* add_instance(ygui::scene* scene);
ygui::object*   add_object(ygui::scene* scene);
ygui::light*    add_light(ygui::scene* scene);

// camera properties
void set_frame(ygui::camera* camera, const frame3f& frame);
void set_lens(ygui::camera* camera, float lens, float aspect, float film);
void set_nearfar(ygui::camera* camera, float near, float far);

// texture properties
void set_texture(
    ygui::texture* texture, const yimg::image<vec4b>& img, bool as_srgb = true);
void set_texture(ygui::texture* texture, const yimg::image<vec4f>& img,
    bool as_float = false);
void set_texture(
    ygui::texture* texture, const yimg::image<vec3b>& img, bool as_srgb = true);
void set_texture(ygui::texture* texture, const yimg::image<vec3f>& img,
    bool as_float = false);
void set_texture(
    ygui::texture* texture, const yimg::image<byte>& img, bool as_srgb = true);
void set_texture(ygui::texture* texture, const yimg::image<float>& img,
    bool as_float = false);

// material properties
void set_emission(ygui::material* material, const vec3f& emission,
    ygui::texture* emission_tex = nullptr);
void set_color(ygui::material* material, const vec3f& color,
    ygui::texture* color_tex = nullptr);
void set_metallic(ygui::material* material, float metallic,
    ygui::texture* metallic_tex = nullptr);
void set_roughness(ygui::material* material, float roughness,
    ygui::texture* roughness_tex = nullptr);
void set_specular(ygui::material* material, float specular,
    ygui::texture* specular_tex = nullptr);
void set_opacity(ygui::material* material, float opacity,
    ygui::texture* opacity_tex = nullptr);
void set_normalmap(ygui::material* material, ygui::texture* normal_tex);

// shape properties
void set_points(ygui::shape* shape, const std::vector<int>& points);
void set_lines(ygui::shape* shape, const std::vector<vec2i>& lines);
void set_triangles(ygui::shape* shape, const std::vector<vec3i>& triangles);
void set_quads(ygui::shape* shape, const std::vector<vec4i>& quads);
void set_positions(ygui::shape* shape, const std::vector<vec3f>& positions);
void set_normals(ygui::shape* shape, const std::vector<vec3f>& normals);
void set_texcoords(ygui::shape* shape, const std::vector<vec2f>& texcoords);
void set_colors(ygui::shape* shape, const std::vector<vec3f>& colors);
void set_tangents(ygui::shape* shape, const std::vector<vec4f>& tangents);

// instance properties
void set_frames(ygui::instance* instance, const std::vector<frame3f>& frames);

// object properties
void set_frame(ygui::object* object, const frame3f& frame);
void set_shape(ygui::object* object, ygui::shape* shape);
void set_material(ygui::object* object, ygui::material* material);
void set_instance(ygui::object* object, ygui::instance* instance);
void set_hidden(ygui::object* object, bool hidden);
void set_highlighted(ygui::object* object, bool highlighted);

// light properties
void set_light(ygui::light* light, const vec3f& position, const vec3f& emission,
    bool directional);

// light size
void clear_lights(ygui::scene* scene);
bool has_max_lights(ygui::scene* scene);

// Draw an OpenGL scene
void draw_scene(ygui::scene* scene, ygui::camera* camera, const vec4i& viewport,
    const scene_params& params);

}  // namespace ygui

// -----------------------------------------------------------------------------
// UI APPLICATION
// -----------------------------------------------------------------------------
namespace ygui {

// Forward declaration of OpenGL window
struct window;

// Input state
struct input {
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
using draw_callback = std::function<void(ygui::window*, const input& input)>;
// Draw callback for drawing widgets
using widgets_callback = std::function<void(ygui::window*, const input& input)>;
// Drop callback that returns that list of dropped strings.
using drop_callback = std::function<void(
    window*, const std::vector<std::string>&, const input& input)>;
// Key callback that returns key codes, pressed/released flag and modifier keys
using key_callback = std::function<void(
    ygui::window*, int key, bool pressed, const input& input)>;
// Char callback that returns ASCII key
using char_callback =
    std::function<void(ygui::window*, unsigned int key, const input& input)>;
// Mouse click callback that returns left/right button, pressed/released flag,
// modifier keys
using click_callback = std::function<void(
    ygui::window*, bool left, bool pressed, const input& input)>;
// Scroll callback that returns scroll amount
using scroll_callback =
    std::function<void(ygui::window*, float amount, const input& input)>;
// Update functions called every frame
using uiupdate_callback =
    std::function<void(ygui::window*, const input& input)>;
// Update functions called every frame
using update_callback = std::function<void(ygui::window*, const input& input)>;

// User interface callcaks
struct ui_callbacks {
  draw_callback     draw_cb     = {};
  widgets_callback  widgets_cb  = {};
  drop_callback     drop_cb     = {};
  key_callback      key_cb      = {};
  char_callback     char_cb     = {};
  click_callback    click_cb    = {};
  scroll_callback   scroll_cb   = {};
  update_callback   update_cb   = {};
  uiupdate_callback uiupdate_cb = {};
};

// run the user interface with the give callbacks
void run_ui(const vec2i& size, const std::string& title,
    const ui_callbacks& callbaks, int widgets_width = 320,
    bool widgets_left = true);

}  // namespace ygui

// -----------------------------------------------------------------------------
// UI WINDOW
// -----------------------------------------------------------------------------
namespace ygui {

// OpenGL window wrapper
struct window {
  GLFWwindow*       win           = nullptr;
  std::string       title         = "";
  draw_callback     draw_cb       = {};
  widgets_callback  widgets_cb    = {};
  drop_callback     drop_cb       = {};
  key_callback      key_cb        = {};
  char_callback     char_cb       = {};
  click_callback    click_cb      = {};
  scroll_callback   scroll_cb     = {};
  update_callback   update_cb     = {};
  uiupdate_callback uiupdate_cb   = {};
  int               widgets_width = 0;
  bool              widgets_left  = true;
  input             input         = {};
  vec4f             background    = {0.15f, 0.15f, 0.15f, 1.0f};
};

// Windows initialization
void init_window(ygui::window* win, const vec2i& size, const std::string& title,
    bool widgets, int widgets_width = 320, bool widgets_left = true);

// Window cleanup
void clear_window(ygui::window* win);

// Set callbacks
void set_draw_callback(ygui::window* win, draw_callback draw_cb);
void set_widgets_callback(ygui::window* win, widgets_callback widgets_cb);
void set_drop_callback(ygui::window* win, drop_callback drop_cb);
void set_key_callback(ygui::window* win, key_callback cb);
void set_char_callback(ygui::window* win, char_callback cb);
void set_click_callback(ygui::window* win, click_callback cb);
void set_scroll_callback(ygui::window* win, scroll_callback cb);
void set_uiupdate_callback(ygui::window* win, uiupdate_callback cb);
void set_update_callback(ygui::window* win, update_callback cb);

// Run loop
void run_ui(ygui::window* win);
void set_close(ygui::window* win, bool close);

}  // namespace ygui

// -----------------------------------------------------------------------------
// OPENGL WIDGETS
// -----------------------------------------------------------------------------
namespace ygui {

bool begin_header(ygui::window* win, const char* title);
void end_header(ygui::window* win);

void draw_label(ygui::window* win, const char* lbl, const std::string& text);

void draw_separator(ygui::window* win);
void continue_line(ygui::window* win);

bool draw_button(ygui::window* win, const char* lbl, bool enabled = true);

bool draw_textinput(ygui::window* win, const char* lbl, std::string& value);

bool draw_slider(
    window* win, const char* lbl, float& value, float min, float max);
bool draw_slider(
    window* win, const char* lbl, vec2f& value, float min, float max);
bool draw_slider(
    window* win, const char* lbl, vec3f& value, float min, float max);
bool draw_slider(
    window* win, const char* lbl, vec4f& value, float min, float max);

bool draw_slider(
    ygui::window* win, const char* lbl, int& value, int min, int max);
bool draw_slider(
    ygui::window* win, const char* lbl, vec2i& value, int min, int max);
bool draw_slider(
    ygui::window* win, const char* lbl, vec3i& value, int min, int max);
bool draw_slider(
    ygui::window* win, const char* lbl, vec4i& value, int min, int max);

bool draw_dragger(ygui::window* win, const char* lbl, float& value,
    float speed = 1.0f, float min = 0.0f, float max = 0.0f);
bool draw_dragger(ygui::window* win, const char* lbl, vec2f& value,
    float speed = 1.0f, float min = 0.0f, float max = 0.0f);
bool draw_dragger(ygui::window* win, const char* lbl, vec3f& value,
    float speed = 1.0f, float min = 0.0f, float max = 0.0f);
bool draw_dragger(ygui::window* win, const char* lbl, vec4f& value,
    float speed = 1.0f, float min = 0.0f, float max = 0.0f);

bool draw_dragger(ygui::window* win, const char* lbl, int& value,
    float speed = 1, int min = 0, int max = 0);
bool draw_dragger(ygui::window* win, const char* lbl, vec2i& value,
    float speed = 1, int min = 0, int max = 0);
bool draw_dragger(ygui::window* win, const char* lbl, vec3i& value,
    float speed = 1, int min = 0, int max = 0);
bool draw_dragger(ygui::window* win, const char* lbl, vec4i& value,
    float speed = 1, int min = 0, int max = 0);

bool draw_checkbox(ygui::window* win, const char* lbl, bool& value);

bool draw_coloredit(ygui::window* win, const char* lbl, vec3f& value);
bool draw_coloredit(ygui::window* win, const char* lbl, vec4f& value);

bool draw_hdrcoloredit(ygui::window* win, const char* lbl, vec3f& value);
bool draw_hdrcoloredit(ygui::window* win, const char* lbl, vec4f& value);

bool draw_combobox(ygui::window* win, const char* lbl, int& idx,
    const std::vector<std::string>& labels);
bool draw_combobox(ygui::window* win, const char* lbl, std::string& value,
    const std::vector<std::string>& labels);
bool draw_combobox(ygui::window* win, const char* lbl, int& idx, int num,
    const std::function<const char*(int)>& labels, bool include_null = false);

template <typename T>
inline bool draw_combobox(ygui::window* win, const char* lbl, int& idx,
    const std::vector<T>& vals, bool include_null = false) {
  return draw_combobox(
      win, lbl, idx, (int)vals.size(),
      [&](int idx) { return vals[idx].name.c_str(); }, include_null);
}
template <typename T>
inline bool draw_combobox(ygui::window* win, const char* lbl, int& idx,
    const std::vector<T*>& vals, bool include_null = false) {
  return draw_combobox(
      win, lbl, idx, (int)vals.size(),
      [&](int idx) { return vals[idx]->name.c_str(); }, include_null);
}
template <typename T>
inline bool draw_combobox(ygui::window* win, const char* lbl, T*& value,
    const std::vector<T*>& vals, bool include_null = false) {
  auto idx = -1;
  for (auto pos = 0; pos < vals.size(); pos++)
    if (vals[pos] == value) idx = pos;
  auto edited = draw_combobox(
      win, lbl, idx, (int)vals.size(),
      [&](int idx) { return vals[idx]->name.c_str(); }, include_null);
  if (edited) {
    value = idx >= 0 ? vals[idx] : nullptr;
  }
  return edited;
}
template <typename T>
inline bool draw_combobox(ygui::window* win, const char* lbl, int& idx,
    const std::vector<std::shared_ptr<T>>& vals, bool include_null = false) {
  return draw_combobox(
      win, lbl, idx, (int)vals.size(),
      [&](int idx) { return vals[idx]->name.c_str(); }, include_null);
}
template <typename T>
inline bool draw_combobox(ygui::window* win, const char* lbl,
    std::shared_ptr<T>& value, const std::vector<std::shared_ptr<T>>& vals,
    bool include_null = false) {
  auto idx = -1;
  for (auto pos = 0; pos < vals.size(); pos++)
    if (vals[pos] == value) idx = pos;
  auto edited = draw_combobox(
      win, lbl, idx, (int)vals.size(),
      [&](int idx) { return vals[idx]->name.c_str(); }, include_null);
  if (edited) {
    value = idx >= 0 ? vals[idx] : nullptr;
  }
  return edited;
}

void draw_progressbar(ygui::window* win, const char* lbl, float fraction);

void draw_histogram(
    ygui::window* win, const char* lbl, const std::vector<float>& values);
void draw_histogram(
    ygui::window* win, const char* lbl, const std::vector<vec2f>& values);
void draw_histogram(
    ygui::window* win, const char* lbl, const std::vector<vec3f>& values);
void draw_histogram(
    ygui::window* win, const char* lbl, const std::vector<vec4f>& values);

bool draw_messages(ygui::window* win);
void push_message(ygui::window* win, const std::string& message);
bool draw_filedialog(ygui::window* win, const char* lbl, std::string& path,
    bool save, const std::string& dirname, const std::string& filename,
    const std::string& filter);
bool draw_filedialog_button(ygui::window* win, const char* button_lbl,
    bool button_active, const char* lbl, std::string& path, bool save,
    const std::string& dirname, const std::string& filename,
    const std::string& filter);

void log_info(ygui::window* win, const std::string& msg);
void log_error(ygui::window* win, const std::string& msg);
void clear_log(ygui::window* win);
void draw_log(ygui::window* win);

}  // namespace ygui

#endif
