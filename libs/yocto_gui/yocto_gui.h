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

#include <yocto/yocto_image.h>
#include <yocto/yocto_math.h>

#include <functional>
#include <memory>
#include <vector>

// forward declaration
struct GLFWwindow;

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::function;
using std::string;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// LOW-LEVEL OPENGL HELPERS
// -----------------------------------------------------------------------------
namespace yocto {

// Commands to setup the opengl context and issue gpu operations.
bool init_ogl(string& error);
void assert_ogl_error();
bool check_ogl_error(string& error);
void clear_ogl_framebuffer(const vec4f& color, bool clear_depth = true);
void set_ogl_viewport(const vec4i& viewport);
void set_ogl_wireframe(bool enabled);
void set_ogl_blending(bool enabled);
void set_ogl_point_size(int size);

// OpenGL texture
struct ogl_texture {
  // Texture properties
  vec2i size      = {0, 0};
  int   nchannels = 0;
  bool  is_srgb   = false;
  bool  is_float  = false;
  bool  linear    = false;
  bool  mipmap    = false;

  // OpenGL state
  uint texture_id = 0;

  // ensuring no copies
  ogl_texture() {}
  ogl_texture(const ogl_texture&) = delete;
  ogl_texture& operator=(ogl_texture&) = delete;
};

// set texture
void set_texture(ogl_texture* texture, const vec2i& size, int nchannels,
    const byte* img, bool as_srgb = false, bool linear = true,
    bool mipmap = true);
void set_texture(ogl_texture* texture, const vec2i& size, int nchannels,
    const float* img, bool as_float = false, bool linear = true,
    bool mipmap = true);

// check if texture is initialized
bool is_initialized(ogl_texture* texture);

// clear texture
void clear_texture(ogl_texture* texture);

// set texture
void set_texture(ogl_texture* texture, const image<vec4b>& img,
    bool as_srgb = true, bool linear = true, bool mipmap = true);
void set_texture(ogl_texture* texture, const image<vec4f>& img,
    bool as_float = false, bool linear = true, bool mipmap = true);
void set_texture(ogl_texture* texture, const image<vec3b>& img,
    bool as_srgb = true, bool linear = true, bool mipmap = true);
void set_texture(ogl_texture* texture, const image<vec3f>& img,
    bool as_float = false, bool linear = true, bool mipmap = true);
void set_texture(ogl_texture* texture, const image<byte>& img,
    bool as_srgb = true, bool linear = true, bool mipmap = true);
void set_texture(ogl_texture* texture, const image<float>& img,
    bool as_float = false, bool linear = true, bool mipmap = true);

// Opengl array/element buffer
struct ogl_arraybuffer {
  // buffer data
  size_t size    = 0;
  int    esize   = 0;
  bool   dynamic = false;
  // OpenGL state
  uint buffer_id = 0;
};

// set buffer
void set_arraybuffer(ogl_arraybuffer* buffer, size_t size, int esize,
    const float* data, bool dynamic = false);

// check if buffer is initialized
bool is_initialized(ogl_arraybuffer* buffer);

// clear buffer
void clear_arraybuffer(ogl_arraybuffer* buffer);

// set buffer
void set_arraybuffer(
    ogl_arraybuffer* buffer, const vector<float>& data, bool dynamic = false);
void set_arraybuffer(
    ogl_arraybuffer* buffer, const vector<vec2f>& data, bool dynamic = false);
void set_arraybuffer(
    ogl_arraybuffer* buffer, const vector<vec3f>& data, bool dynamic = false);
void set_arraybuffer(
    ogl_arraybuffer* buffer, const vector<vec4f>& data, bool dynamic = false);

// Opengl draw elements
enum struct ogl_element_type { points, lines, triangles };

// Opengl array/element buffer
struct ogl_elementbuffer {
  // buffer data
  size_t           size    = 0;
  ogl_element_type element = ogl_element_type::points;
  bool             dynamic = false;
  // OpenGL state
  uint buffer_id = 0;
};

// set buffer
void set_elementbuffer(ogl_elementbuffer* buffer, size_t size,
    ogl_element_type element, const int* data, bool dynamic = false);

// check if buffer is initialized
bool is_initialized(ogl_elementbuffer* buffer);

// clear buffer
void clear_elementbuffer(ogl_elementbuffer* buffer);

// set buffer
void set_elementbuffer(
    ogl_elementbuffer* buffer, const vector<int>& points, bool dynamic = false);
void set_elementbuffer(ogl_elementbuffer* buffer, const vector<vec2i>& lines,
    bool dynamic = false);
void set_elementbuffer(ogl_elementbuffer* buffer,
    const vector<vec3i>& triangles, bool dynamic = false);

// Opengl program
struct ogl_program {
  // program code
  string vertex_code;
  string fragment_code;
  // OpenGL state
  uint program_id  = 0;
  uint vertex_id   = 0;
  uint fragment_id = 0;
  uint array_id    = 0;
};

// initialize program
bool init_program(ogl_program* program, const string& vertex,
    const string& fragment, string& error, string& errorlog);
bool is_initialized(const ogl_program* program);

// clear program
void clear_program(ogl_program* program);

// bind program
void bind_program(ogl_program* program);
// unbind program
void unbind_program(ogl_program* program);
// unbind program
void unbind_program();

// get uniform location
int get_uniform_location(ogl_program* program, const char* name);

// set uniforms
void set_uniform(ogl_program* program, int location, int value);
void set_uniform(ogl_program* program, int location, const vec2i& value);
void set_uniform(ogl_program* program, int location, const vec3i& value);
void set_uniform(ogl_program* program, int location, const vec4i& value);
void set_uniform(ogl_program* program, int location, float value);
void set_uniform(ogl_program* program, int location, const vec2f& value);
void set_uniform(ogl_program* program, int location, const vec3f& value);
void set_uniform(ogl_program* program, int location, const vec4f& value);
void set_uniform(ogl_program* program, int location, const mat2f& value);
void set_uniform(ogl_program* program, int location, const mat3f& value);
void set_uniform(ogl_program* program, int location, const mat4f& value);
void set_uniform(ogl_program* program, int location, const frame2f& value);
void set_uniform(ogl_program* program, int location, const frame3f& value);
template <typename T>
inline void set_uniform(
    ogl_program* program, const char* name, const T& value) {
  return set_uniform(program, get_uniform_location(program, name), value);
}

// set uniform texture
void set_uniform(
    ogl_program* program, int location, const ogl_texture* texture, int unit);
void set_uniform(ogl_program* program, const char* name,
    const ogl_texture* texture, int unit);
void set_uniform(ogl_program* program, int location, int location_on,
    const ogl_texture* texture, int unit);
void set_uniform(ogl_program* program, const char* name, const char* name_on,
    const ogl_texture* texture, int unit);

// get attribute location
int get_attribute_location(ogl_program* program, const char* name);

// set vertex attributes
void set_attribute(ogl_program* program, int location, ogl_arraybuffer* buffer);
void set_attribute(
    ogl_program* program, const char* name, ogl_arraybuffer* buffer);

// set vertex attributes
void set_attribute(ogl_program* program, int location, float value);
void set_attribute(ogl_program* program, int location, const vec2f& value);
void set_attribute(ogl_program* program, int location, const vec3f& value);
void set_attribute(ogl_program* program, int location, const vec4f& value);
template <typename T>
inline void set_attribute(
    ogl_program* program, const char* name, const T& value) {
  return set_attribute(program, get_attribute_location(program, name), value);
}

// set vertex attributes
template <typename T>
inline void set_attribute(ogl_program* program, int location,
    ogl_arraybuffer* buffer, const T& value) {
  if (buffer && is_initialized(buffer)) {
    return set_attribute(program, location, buffer);
  } else {
    set_attribute(program, location, value);
  }
}
template <typename T>
inline void set_attribute(ogl_program* program, const char* name,
    ogl_arraybuffer* buffer, const T& def) {
  set_attribute(program, get_attribute_location(program, name), buffer, def);
}

// draw elements
void draw_elements(ogl_elementbuffer* buffer);

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE DRAWING
// -----------------------------------------------------------------------------
namespace yocto {

// OpenGL image data
struct ogl_image {
  ogl_image() {}
  ogl_image(const ogl_image&) = delete;
  ogl_image& operator=(const ogl_image&) = delete;
  ~ogl_image();

  ogl_program*       program   = new ogl_program{};
  ogl_texture*       texture   = new ogl_texture{};
  ogl_arraybuffer*   texcoords = new ogl_arraybuffer{};
  ogl_elementbuffer* triangles = new ogl_elementbuffer{};
};

// create image drawing program
bool init_image(ogl_image* oimg);
bool is_initialized(const ogl_image* oimg);

// clear image
void clear_image(ogl_image* oimg);

// update image data
void set_image(ogl_image* oimg, const image<vec4f>& img, bool linear = false,
    bool mipmap = false);
void set_image(ogl_image* oimg, const image<vec4b>& img, bool linear = false,
    bool mipmap = false);

// OpenGL image drawing params
struct ogl_image_params {
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
void draw_image(ogl_image* image, const ogl_image_params& params);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE DRAWING
// -----------------------------------------------------------------------------
namespace yocto {

// Opengl caemra
struct ogl_camera {
  frame3f frame    = identity3x4f;
  float   lens     = 0.050;
  float   aspect   = 1.000;
  float   film     = 0.036;
  float   near     = 0.001;
  float   far      = 10000;
  float   aperture = 0;
  float   focus    = 0;
};

// Opengl material
struct ogl_material {
  // material
  vec3f        emission      = {0, 0, 0};
  vec3f        color         = {0, 0, 0};
  float        metallic      = 0;
  float        roughness     = 0;
  float        specular      = 0;
  float        opacity       = 1;
  ogl_texture* emission_tex  = nullptr;
  ogl_texture* color_tex     = nullptr;
  ogl_texture* metallic_tex  = nullptr;
  ogl_texture* roughness_tex = nullptr;
  ogl_texture* specular_tex  = nullptr;
  ogl_texture* opacity_tex   = nullptr;
  ogl_texture* normal_tex    = nullptr;
};

// Opengl shape
struct ogl_shape {
  // vertex buffers
  ogl_arraybuffer*   positions      = new ogl_arraybuffer{};
  ogl_arraybuffer*   normals        = new ogl_arraybuffer{};
  ogl_arraybuffer*   texcoords      = new ogl_arraybuffer{};
  ogl_arraybuffer*   colors         = new ogl_arraybuffer{};
  ogl_arraybuffer*   tangents       = new ogl_arraybuffer{};
  ogl_elementbuffer* points         = new ogl_elementbuffer{};
  ogl_elementbuffer* lines          = new ogl_elementbuffer{};
  ogl_elementbuffer* triangles      = new ogl_elementbuffer{};
  ogl_elementbuffer* quads          = new ogl_elementbuffer{};
  ogl_elementbuffer* edges          = new ogl_elementbuffer{};
  float              points_size    = 10;
  float              line_thickness = 4;

  ogl_shape() {}
  ogl_shape(const ogl_shape&) = delete;
  ogl_shape& operator=(const ogl_shape&) = delete;
  ~ogl_shape();
};

// Opengl instance
struct ogl_instance {
  // instance properties
  frame3f       frame       = identity3x4f;
  ogl_shape*    shape       = nullptr;
  ogl_material* material    = nullptr;
  bool          hidden      = false;
  bool          highlighted = false;
};

// Light type
enum struct ogl_light_type { point = 0, directional };

// Opengl light
struct ogl_light {
  vec3f          position = {0, 0, 0};
  vec3f          emission = {0, 0, 0};
  ogl_light_type type     = ogl_light_type::point;
  bool           camera   = false;
};

// Opengl scene
struct ogl_scene {
  ogl_scene() {}
  ogl_scene(const ogl_scene&) = delete;
  ogl_scene& operator=(const ogl_scene&) = delete;
  ~ogl_scene();

  // scene objects
  vector<ogl_camera*>   cameras   = {};
  vector<ogl_instance*> instances = {};
  vector<ogl_shape*>    shapes    = {};
  vector<ogl_material*> materials = {};
  vector<ogl_texture*>  textures  = {};
  vector<ogl_light*>    lights    = {};

  // OpenGL state
  ogl_program* program = new ogl_program{};
};

// Shading type
enum struct ogl_shading_type { lights, eyelight, camlights };

// Shading name
const auto ogl_shading_names = vector<string>{
    "lights", "eyelight", "camlights"};

// Draw options
struct ogl_scene_params {
  int              resolution       = 1280;
  bool             wireframe        = false;
  bool             edges            = false;
  float            edge_offset      = 0.01f;
  ogl_shading_type shading          = ogl_shading_type::camlights;
  float            exposure         = 0;
  float            gamma            = 2.2f;
  vec3f            ambient          = {0, 0, 0};
  bool             double_sided     = true;
  bool             non_rigid_frames = true;
  float            near             = 0.01f;
  float            far              = 10000.0f;
  vec4f            background       = vec4f{0.15f, 0.15f, 0.15f, 1.0f};
};

// Initialize an OpenGL scene
void init_scene(ogl_scene* scene);
bool is_initialized(const ogl_scene* scene);

// Clear an OpenGL scene
void clear_scene(ogl_scene* scene);

// add scene elements
ogl_camera*   add_camera(ogl_scene* scene);
ogl_texture*  add_texture(ogl_scene* scene);
ogl_material* add_material(ogl_scene* scene);
ogl_shape*    add_shape(ogl_scene* scene);
ogl_instance* add_instance(ogl_scene* scene);
ogl_light*    add_light(ogl_scene* scene);

// camera properties
void set_frame(ogl_camera* camera, const frame3f& frame);
void set_lens(ogl_camera* camera, float lens, float aspect, float film);
void set_nearfar(ogl_camera* camera, float near, float far);

// material properties
void set_emission(ogl_material* material, const vec3f& emission,
    ogl_texture* emission_tex = nullptr);
void set_color(ogl_material* material, const vec3f& color,
    ogl_texture* color_tex = nullptr);
void set_metallic(ogl_material* material, float metallic,
    ogl_texture* metallic_tex = nullptr);
void set_roughness(ogl_material* material, float roughness,
    ogl_texture* roughness_tex = nullptr);
void set_specular(ogl_material* material, float specular,
    ogl_texture* specular_tex = nullptr);
void set_opacity(
    ogl_material* material, float opacity, ogl_texture* opacity_tex = nullptr);
void set_normalmap(ogl_material* material, ogl_texture* normal_tex);

// shape properties
void set_points(ogl_shape* shape, const vector<int>& points);
void set_lines(ogl_shape* shape, const vector<vec2i>& lines);
void set_triangles(ogl_shape* shape, const vector<vec3i>& triangles);
void set_quads(ogl_shape* shape, const vector<vec4i>& quads);
void set_edges(ogl_shape* shape, const vector<vec3i>& triangles,
    const vector<vec4i>& quads);
void set_positions(ogl_shape* shape, const vector<vec3f>& positions);
void set_normals(ogl_shape* shape, const vector<vec3f>& normals);
void set_texcoords(ogl_shape* shape, const vector<vec2f>& texcoords);
void set_colors(ogl_shape* shape, const vector<vec3f>& colors);
void set_tangents(ogl_shape* shape, const vector<vec4f>& tangents);

// instance properties
void set_frame(ogl_instance* instance, const frame3f& frame);
void set_shape(ogl_instance* instance, ogl_shape* shape);
void set_material(ogl_instance* instance, ogl_material* material);
void set_hidden(ogl_instance* instance, bool hidden);
void set_highlighted(ogl_instance* instance, bool highlighted);

// shortcuts
ogl_camera*   add_camera(ogl_scene* scene, const frame3f& frame, float lens,
      float aspect, float film = 0.036, float near = 0.001, float far = 10000);
ogl_material* add_material(ogl_scene* scene, const vec3f& emission,
    const vec3f& color, float specular, float metallic, float roughness,
    ogl_texture* emission_tex = nullptr, ogl_texture* color_tex = nullptr,
    ogl_texture* specular_tex = nullptr, ogl_texture* metallic_tex = nullptr,
    ogl_texture* roughness_tex = nullptr, ogl_texture* normalmap_tex = nullptr);
ogl_shape*    add_shape(ogl_scene* scene, const vector<int>& points,
       const vector<vec2i>& lines, const vector<vec3i>& triangles,
       const vector<vec4i>& quads, const vector<vec3f>& positions,
       const vector<vec3f>& normals, const vector<vec2f>& texcoords,
       const vector<vec3f>& colors, bool edges = false);
ogl_instance* add_instance(ogl_scene* scene, const frame3f& frame,
    ogl_shape* shape, ogl_material* material, bool hidden = false,
    bool highlighted = false);

// light properties
void add_default_lights(ogl_scene* scene);
void set_light(ogl_light* light, const vec3f& position, const vec3f& emission,
    ogl_light_type type, bool camera);

// light size
void clear_lights(ogl_scene* scene);
bool has_max_lights(ogl_scene* scene);

// Draw an OpenGL scene
void draw_scene(ogl_scene* scene, ogl_camera* camera, const vec4i& viewport,
    const ogl_scene_params& params);

}  // namespace yocto

// -----------------------------------------------------------------------------
// UI APPLICATION
// -----------------------------------------------------------------------------
namespace yocto {

// Forward declaration of OpenGL window
struct gui_window;

// Input state
struct gui_input {
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

// Init callback called after the window has opened
using init_callback = function<void(gui_window*, const gui_input& input)>;
// Clear callback called after the window is cloased
using clear_callback = function<void(gui_window*, const gui_input& input)>;
// Draw callback called every frame and when resizing
using draw_callback = function<void(gui_window*, const gui_input& input)>;
// Draw callback for drawing widgets
using widgets_callback = function<void(gui_window*, const gui_input& input)>;
// Drop callback that returns that list of dropped strings.
using drop_callback =
    function<void(gui_window*, const vector<string>&, const gui_input& input)>;
// Key callback that returns key codes, pressed/released flag and modifier keys
using key_callback =
    function<void(gui_window*, int key, bool pressed, const gui_input& input)>;
// Char callback that returns ASCII key
using char_callback =
    function<void(gui_window*, unsigned int key, const gui_input& input)>;
// Mouse click callback that returns left/right button, pressed/released flag,
// modifier keys
using click_callback = function<void(
    gui_window*, bool left, bool pressed, const gui_input& input)>;
// Scroll callback that returns scroll amount
using scroll_callback =
    function<void(gui_window*, float amount, const gui_input& input)>;
// Update functions called every frame
using uiupdate_callback = function<void(gui_window*, const gui_input& input)>;
// Update functions called every frame
using update_callback = function<void(gui_window*, const gui_input& input)>;

// User interface callcaks
struct gui_callbacks {
  init_callback     init_cb     = {};
  clear_callback    clear_cb    = {};
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
void run_ui(const vec2i& size, const string& title,
    const gui_callbacks& callbaks, int widgets_width = 320,
    bool widgets_left = true);

}  // namespace yocto

// -----------------------------------------------------------------------------
// UI WINDOW
// -----------------------------------------------------------------------------
namespace yocto {

// OpenGL window wrapper
struct gui_window {
  GLFWwindow*       win           = nullptr;
  string            title         = "";
  init_callback     init_cb       = {};
  clear_callback    clear_cb      = {};
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
  gui_input         input         = {};
  vec4f             background    = {0.15f, 0.15f, 0.15f, 1.0f};
};

// Windows initialization
void init_window(gui_window* win, const vec2i& size, const string& title,
    bool widgets, int widgets_width = 320, bool widgets_left = true);

// Window cleanup
void clear_window(gui_window* win);

// Set callbacks
void set_init_callback(gui_window* win, init_callback init_cb);
void set_clear_callback(gui_window* win, clear_callback clear_cb);
void set_draw_callback(gui_window* win, draw_callback draw_cb);
void set_widgets_callback(gui_window* win, widgets_callback widgets_cb);
void set_drop_callback(gui_window* win, drop_callback drop_cb);
void set_key_callback(gui_window* win, key_callback cb);
void set_char_callback(gui_window* win, char_callback cb);
void set_click_callback(gui_window* win, click_callback cb);
void set_scroll_callback(gui_window* win, scroll_callback cb);
void set_uiupdate_callback(gui_window* win, uiupdate_callback cb);
void set_update_callback(gui_window* win, update_callback cb);

// Run loop
void run_ui(gui_window* win);
void set_close(gui_window* win, bool close);

}  // namespace yocto

// -----------------------------------------------------------------------------
// OPENGL WIDGETS
// -----------------------------------------------------------------------------
namespace yocto {

bool begin_header(gui_window* win, const char* title);
void end_header(gui_window* win);

void draw_label(gui_window* win, const char* lbl, const string& text);

void draw_separator(gui_window* win);
void continue_line(gui_window* win);

bool draw_button(gui_window* win, const char* lbl, bool enabled = true);

bool draw_textinput(gui_window* win, const char* lbl, string& value);

bool draw_slider(
    gui_window* win, const char* lbl, float& value, float min, float max);
bool draw_slider(
    gui_window* win, const char* lbl, vec2f& value, float min, float max);
bool draw_slider(
    gui_window* win, const char* lbl, vec3f& value, float min, float max);
bool draw_slider(
    gui_window* win, const char* lbl, vec4f& value, float min, float max);

bool draw_slider(
    gui_window* win, const char* lbl, int& value, int min, int max);
bool draw_slider(
    gui_window* win, const char* lbl, vec2i& value, int min, int max);
bool draw_slider(
    gui_window* win, const char* lbl, vec3i& value, int min, int max);
bool draw_slider(
    gui_window* win, const char* lbl, vec4i& value, int min, int max);

bool draw_dragger(gui_window* win, const char* lbl, float& value,
    float speed = 1.0f, float min = 0.0f, float max = 0.0f);
bool draw_dragger(gui_window* win, const char* lbl, vec2f& value,
    float speed = 1.0f, float min = 0.0f, float max = 0.0f);
bool draw_dragger(gui_window* win, const char* lbl, vec3f& value,
    float speed = 1.0f, float min = 0.0f, float max = 0.0f);
bool draw_dragger(gui_window* win, const char* lbl, vec4f& value,
    float speed = 1.0f, float min = 0.0f, float max = 0.0f);

bool draw_dragger(gui_window* win, const char* lbl, int& value, float speed = 1,
    int min = 0, int max = 0);
bool draw_dragger(gui_window* win, const char* lbl, vec2i& value,
    float speed = 1, int min = 0, int max = 0);
bool draw_dragger(gui_window* win, const char* lbl, vec3i& value,
    float speed = 1, int min = 0, int max = 0);
bool draw_dragger(gui_window* win, const char* lbl, vec4i& value,
    float speed = 1, int min = 0, int max = 0);

bool draw_checkbox(gui_window* win, const char* lbl, bool& value);
bool draw_checkbox(gui_window* win, const char* lbl, bool& value, bool invert);

bool draw_coloredit(gui_window* win, const char* lbl, vec3f& value);
bool draw_coloredit(gui_window* win, const char* lbl, vec4f& value);

bool draw_hdrcoloredit(gui_window* win, const char* lbl, vec3f& value);
bool draw_hdrcoloredit(gui_window* win, const char* lbl, vec4f& value);

bool draw_combobox(
    gui_window* win, const char* lbl, int& idx, const vector<string>& labels);
bool draw_combobox(gui_window* win, const char* lbl, string& value,
    const vector<string>& labels);
bool draw_combobox(gui_window* win, const char* lbl, int& idx, int num,
    const function<string(int)>& labels, bool include_null = false);

template <typename T>
inline bool draw_combobox(gui_window* win, const char* lbl, T*& value,
    const vector<T*>& vals, bool include_null = false) {
  auto idx = -1;
  for (auto pos = 0; pos < vals.size(); pos++)
    if (vals[pos] == value) idx = pos;
  auto edited = draw_combobox(
      win, lbl, idx, (int)vals.size(), [&](int idx) { return vals[idx]->name; },
      include_null);
  if (edited) {
    value = idx >= 0 ? vals[idx] : nullptr;
  }
  return edited;
}

template <typename T>
inline bool draw_combobox(gui_window* win, const char* lbl, T*& value,
    const vector<T*>& vals, const vector<string>& labels,
    bool include_null = false) {
  auto idx = -1;
  for (auto pos = 0; pos < vals.size(); pos++)
    if (vals[pos] == value) idx = pos;
  auto edited = draw_combobox(
      win, lbl, idx, (int)vals.size(), [&](int idx) { return labels[idx]; },
      include_null);
  if (edited) {
    value = idx >= 0 ? vals[idx] : nullptr;
  }
  return edited;
}

void draw_progressbar(gui_window* win, const char* lbl, float fraction);
void draw_progressbar(gui_window* win, const char* lbl, int current, int total);

void draw_histogram(
    gui_window* win, const char* lbl, const vector<float>& values);
void draw_histogram(
    gui_window* win, const char* lbl, const vector<vec2f>& values);
void draw_histogram(
    gui_window* win, const char* lbl, const vector<vec3f>& values);
void draw_histogram(
    gui_window* win, const char* lbl, const vector<vec4f>& values);

bool draw_messages(gui_window* win);
void push_message(gui_window* win, const string& message);
bool draw_filedialog(gui_window* win, const char* lbl, string& path, bool save,
    const string& dirname, const string& filename, const string& filter);
bool draw_filedialog_button(gui_window* win, const char* button_lbl,
    bool button_active, const char* lbl, string& path, bool save,
    const string& dirname, const string& filename, const string& filter);

void log_info(gui_window* win, const string& msg);
void log_error(gui_window* win, const string& msg);
void clear_log(gui_window* win);
void draw_log(gui_window* win);

}  // namespace yocto

#endif
