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
// ALIASES
// -----------------------------------------------------------------------------
namespace yocto::gui {

// Namespace aliases
namespace gui = yocto::gui;
namespace img = yocto::image;

}  // namespace yocto::gui

// -----------------------------------------------------------------------------
// LOW-LEVEL OPENGL HELPERS
// -----------------------------------------------------------------------------
namespace yocto::gui {

// Commands to setup the opengl context and issue gpu operations.
bool init_opengl(std::string& error);
void assert_error();
bool check_error(std::string& error);
void clear_framebuffer(const vec4f& color, bool clear_depth = true);
void set_viewport(const vec4i& viewport);
void set_wireframe(bool enabled);
void set_blending(bool enabled);
void set_point_size(int size);

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
void set_texture(gui::ogl_texture* texture, const vec2i& size, int nchannels,
    const byte* img, bool as_srgb = false, bool linear = true,
    bool mipmap = true);
void set_texture(gui::ogl_texture* texture, const vec2i& size, int nchannels,
    const float* img, bool as_float = false, bool linear = true,
    bool mipmap = true);

// check if texture is initialized
bool is_initialized(gui::ogl_texture* texture);

// clear texture
void clear_texture(gui::ogl_texture* texture);

// set texture
void set_texture(gui::ogl_texture* texture, const img::image<vec4b>& img,
    bool as_srgb = true, bool linear = true, bool mipmap = true);
void set_texture(gui::ogl_texture* texture, const img::image<vec4f>& img,
    bool as_float = false, bool linear = true, bool mipmap = true);
void set_texture(gui::ogl_texture* texture, const img::image<vec3b>& img,
    bool as_srgb = true, bool linear = true, bool mipmap = true);
void set_texture(gui::ogl_texture* texture, const img::image<vec3f>& img,
    bool as_float = false, bool linear = true, bool mipmap = true);
void set_texture(gui::ogl_texture* texture, const img::image<byte>& img,
    bool as_srgb = true, bool linear = true, bool mipmap = true);
void set_texture(gui::ogl_texture* texture, const img::image<float>& img,
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
void set_arraybuffer(gui::ogl_arraybuffer* buffer, size_t size, int esize,
    const float* data, bool dynamic = false);

// check if buffer is initialized
bool is_initialized(gui::ogl_arraybuffer* buffer);

// clear buffer
void clear_arraybuffer(gui::ogl_arraybuffer* buffer);

// set buffer
void set_arraybuffer(gui::ogl_arraybuffer* buffer,
    const std::vector<float>& data, bool dynamic = false);
void set_arraybuffer(gui::ogl_arraybuffer* buffer,
    const std::vector<vec2f>& data, bool dynamic = false);
void set_arraybuffer(gui::ogl_arraybuffer* buffer,
    const std::vector<vec3f>& data, bool dynamic = false);
void set_arraybuffer(gui::ogl_arraybuffer* buffer,
    const std::vector<vec4f>& data, bool dynamic = false);

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
void set_elementbuffer(gui::ogl_elementbuffer* buffer, size_t size,
    ogl_element_type element, const int* data, bool dynamic = false);

// check if buffer is initialized
bool is_initialized(gui::ogl_elementbuffer* buffer);

// clear buffer
void clear_elementbuffer(gui::ogl_elementbuffer* buffer);

// set buffer
void set_elementbuffer(gui::ogl_elementbuffer* buffer,
    const std::vector<int>& points, bool dynamic = false);
void set_elementbuffer(gui::ogl_elementbuffer* buffer,
    const std::vector<vec2i>& lines, bool dynamic = false);
void set_elementbuffer(gui::ogl_elementbuffer* buffer,
    const std::vector<vec3i>& triangles, bool dynamic = false);

// Opengl program
struct ogl_program {
  // program code
  std::string vertex_code;
  std::string fragment_code;
  // OpenGL state
  uint program_id  = 0;
  uint vertex_id   = 0;
  uint fragment_id = 0;
  uint array_id    = 0;
};

// initialize program
bool init_program(gui::ogl_program* program, const std::string& vertex,
    const std::string& fragment, std::string& error, std::string& errorlog);
bool is_initialized(const gui::ogl_program* program);

// clear program
void clear_program(gui::ogl_program* program);

// bind program
void bind_program(gui::ogl_program* program);
// unbind program
void unbind_program(gui::ogl_program* program);
// unbind program
void unbind_program();

// get uniform location
int get_uniform_location(gui::ogl_program* program, const char* name);

// set uniforms
void set_uniform(gui::ogl_program* program, int location, int value);
void set_uniform(gui::ogl_program* program, int location, const vec2i& value);
void set_uniform(gui::ogl_program* program, int location, const vec3i& value);
void set_uniform(gui::ogl_program* program, int location, const vec4i& value);
void set_uniform(gui::ogl_program* program, int location, float value);
void set_uniform(gui::ogl_program* program, int location, const vec2f& value);
void set_uniform(gui::ogl_program* program, int location, const vec3f& value);
void set_uniform(gui::ogl_program* program, int location, const vec4f& value);
void set_uniform(gui::ogl_program* program, int location, const mat2f& value);
void set_uniform(gui::ogl_program* program, int location, const mat3f& value);
void set_uniform(gui::ogl_program* program, int location, const mat4f& value);
void set_uniform(gui::ogl_program* program, int location, const frame2f& value);
void set_uniform(gui::ogl_program* program, int location, const frame3f& value);
template <typename T>
inline void set_uniform(
    gui::ogl_program* program, const char* name, const T& value) {
  return set_uniform(program, get_uniform_location(program, name), value);
}

// set uniform texture
void set_uniform(gui::ogl_program* program, int location,
    const gui::ogl_texture* texture, int unit);
void set_uniform(gui::ogl_program* program, const char* name,
    const gui::ogl_texture* texture, int unit);
void set_uniform(gui::ogl_program* program, int location, int location_on,
    const gui::ogl_texture* texture, int unit);
void set_uniform(gui::ogl_program* program, const char* name,
    const char* name_on, const gui::ogl_texture* texture, int unit);

// get attribute location
int get_attribute_location(gui::ogl_program* program, const char* name);

// set vertex attributes
void set_attribute(
    gui::ogl_program* program, int location, gui::ogl_arraybuffer* buffer);
void set_attribute(
    gui::ogl_program* program, const char* name, gui::ogl_arraybuffer* buffer);

// set vertex attributes
void set_attribute(gui::ogl_program* program, int location, float value);
void set_attribute(gui::ogl_program* program, int location, const vec2f& value);
void set_attribute(gui::ogl_program* program, int location, const vec3f& value);
void set_attribute(gui::ogl_program* program, int location, const vec4f& value);
template <typename T>
inline void set_attribute(
    gui::ogl_program* program, const char* name, const T& value) {
  return set_attribute(program, get_attribute_location(program, name), value);
}

// set vertex attributes
template <typename T>
inline void set_attribute(gui::ogl_program* program, int location,
    gui::ogl_arraybuffer* buffer, const T& value) {
  if (buffer && is_initialized(buffer)) {
    return set_attribute(program, location, buffer);
  } else {
    set_attribute(program, location, value);
  }
}
template <typename T>
inline void set_attribute(gui::ogl_program* program, const char* name,
    gui::ogl_arraybuffer* buffer, const T& def) {
  set_attribute(program, get_attribute_location(program, name), buffer, def);
}

// draw elements
void draw_elements(gui::ogl_elementbuffer* buffer);

}  // namespace yocto::gui

// -----------------------------------------------------------------------------
// IMAGE DRAWING
// -----------------------------------------------------------------------------
namespace yocto::gui {

// OpenGL image data
struct ogl_image {
  ogl_image() {}
  ogl_image(const ogl_image&) = delete;
  ogl_image& operator=(const ogl_image&) = delete;
  ~ogl_image();

  gui::ogl_program*       program   = new gui::ogl_program{};
  gui::ogl_texture*       texture   = new gui::ogl_texture{};
  gui::ogl_arraybuffer*   texcoords = new gui::ogl_arraybuffer{};
  gui::ogl_elementbuffer* triangles = new gui::ogl_elementbuffer{};
};

// create image drawing program
bool init_image(gui::ogl_image* image);
bool is_initialized(const gui::ogl_image* image);

// clear image
void clear_image(gui::ogl_image* image);

// update image data
void set_image(gui::ogl_image* image, const img::image<vec4f>& img,
    bool linear = false, bool mipmap = false);
void set_image(gui::ogl_image* image, const img::image<vec4b>& img,
    bool linear = false, bool mipmap = false);

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
void draw_image(gui::ogl_image* image, const ogl_image_params& params);

}  // namespace yocto::gui

// -----------------------------------------------------------------------------
// SCENE DRAWING
// -----------------------------------------------------------------------------
namespace yocto::gui {

// Opengl caemra
struct ogl_camera {
  frame3f frame  = identity3x4f;
  float   lens   = 0.050;
  float   aspect = 1.000;
  float   film   = 0.036;
  float   near   = 0.001;
  float   far    = 10000;
};

// Opengl material
struct ogl_material {
  // material
  vec3f             emission      = {0, 0, 0};
  vec3f             color         = {0, 0, 0};
  float             metallic      = 0;
  float             roughness     = 0;
  float             specular      = 0;
  float             opacity       = 1;
  gui::ogl_texture* emission_tex  = nullptr;
  gui::ogl_texture* color_tex     = nullptr;
  gui::ogl_texture* metallic_tex  = nullptr;
  gui::ogl_texture* roughness_tex = nullptr;
  gui::ogl_texture* specular_tex  = nullptr;
  gui::ogl_texture* opacity_tex   = nullptr;
  gui::ogl_texture* normal_tex    = nullptr;
};

// Opengl shape
struct ogl_shape {
  // vertex buffers
  gui::ogl_arraybuffer*   positions      = new gui::ogl_arraybuffer{};
  gui::ogl_arraybuffer*   normals        = new gui::ogl_arraybuffer{};
  gui::ogl_arraybuffer*   texcoords      = new gui::ogl_arraybuffer{};
  gui::ogl_arraybuffer*   colors         = new gui::ogl_arraybuffer{};
  gui::ogl_arraybuffer*   tangents       = new gui::ogl_arraybuffer{};
  gui::ogl_elementbuffer* points         = new gui::ogl_elementbuffer{};
  gui::ogl_elementbuffer* lines          = new gui::ogl_elementbuffer{};
  gui::ogl_elementbuffer* triangles      = new gui::ogl_elementbuffer{};
  gui::ogl_elementbuffer* quads          = new gui::ogl_elementbuffer{};
  gui::ogl_elementbuffer* edges          = new gui::ogl_elementbuffer{};
  float                   points_size    = 10;
  float                   line_thickness = 4;

  ogl_shape() {}
  ogl_shape(const ogl_shape&) = delete;
  ogl_shape& operator=(const ogl_shape&) = delete;
  ~ogl_shape();
};

// Opengl instance
struct ogl_instance {
  std::vector<frame3f> frames = {};
};

// Opengl object
struct ogl_object {
  // object properties
  frame3f            frame       = identity3x4f;
  gui::ogl_shape*    shape       = nullptr;
  gui::ogl_material* material    = nullptr;
  gui::ogl_instance* instance    = nullptr;
  bool               hidden      = false;
  bool               highlighted = false;
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
  std::vector<gui::ogl_camera*>   cameras   = {};
  std::vector<gui::ogl_object*>   objects   = {};
  std::vector<gui::ogl_shape*>    shapes    = {};
  std::vector<gui::ogl_material*> materials = {};
  std::vector<gui::ogl_instance*> instances = {};
  std::vector<gui::ogl_texture*>  textures  = {};
  std::vector<gui::ogl_light*>    lights    = {};

  // OpenGL state
  gui::ogl_program* program = new gui::ogl_program{};
};

// Shading type
enum struct ogl_shading_type { lights, eyelight, camlights };

// Shading name
const auto ogl_shading_names = std::vector<std::string>{
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
void init_scene(gui::ogl_scene* scene);
bool is_initialized(const gui::ogl_scene* scene);

// Clear an OpenGL scene
void clear_scene(gui::ogl_scene* scene);

// add scene elements
gui::ogl_camera*   add_camera(gui::ogl_scene* scene);
gui::ogl_texture*  add_texture(gui::ogl_scene* scene);
gui::ogl_material* add_material(gui::ogl_scene* scene);
gui::ogl_shape*    add_shape(gui::ogl_scene* scene);
gui::ogl_instance* add_instance(gui::ogl_scene* scene);
gui::ogl_object*   add_object(gui::ogl_scene* scene);
gui::ogl_light*    add_light(gui::ogl_scene* scene);

// camera properties
void set_frame(gui::ogl_camera* camera, const frame3f& frame);
void set_lens(gui::ogl_camera* camera, float lens, float aspect, float film);
void set_nearfar(gui::ogl_camera* camera, float near, float far);

// material properties
void set_emission(gui::ogl_material* material, const vec3f& emission,
    gui::ogl_texture* emission_tex = nullptr);
void set_color(gui::ogl_material* material, const vec3f& color,
    gui::ogl_texture* color_tex = nullptr);
void set_metallic(gui::ogl_material* material, float metallic,
    gui::ogl_texture* metallic_tex = nullptr);
void set_roughness(gui::ogl_material* material, float roughness,
    gui::ogl_texture* roughness_tex = nullptr);
void set_specular(gui::ogl_material* material, float specular,
    gui::ogl_texture* specular_tex = nullptr);
void set_opacity(gui::ogl_material* material, float opacity,
    gui::ogl_texture* opacity_tex = nullptr);
void set_normalmap(gui::ogl_material* material, gui::ogl_texture* normal_tex);

// shape properties
void set_points(gui::ogl_shape* shape, const std::vector<int>& points);
void set_lines(gui::ogl_shape* shape, const std::vector<vec2i>& lines);
void set_triangles(gui::ogl_shape* shape, const std::vector<vec3i>& triangles);
void set_quads(gui::ogl_shape* shape, const std::vector<vec4i>& quads);
void set_edges(gui::ogl_shape* shape, const std::vector<vec3i>& triangles,
    const std::vector<vec4i>& quads);
void set_positions(gui::ogl_shape* shape, const std::vector<vec3f>& positions);
void set_normals(gui::ogl_shape* shape, const std::vector<vec3f>& normals);
void set_texcoords(gui::ogl_shape* shape, const std::vector<vec2f>& texcoords);
void set_colors(gui::ogl_shape* shape, const std::vector<vec3f>& colors);
void set_tangents(gui::ogl_shape* shape, const std::vector<vec4f>& tangents);

// instance properties
void set_frames(
    gui::ogl_instance* instance, const std::vector<frame3f>& frames);

// object properties
void set_frame(gui::ogl_object* object, const frame3f& frame);
void set_shape(gui::ogl_object* object, gui::ogl_shape* shape);
void set_material(gui::ogl_object* object, gui::ogl_material* material);
void set_instance(gui::ogl_object* object, gui::ogl_instance* instance);
void set_hidden(gui::ogl_object* object, bool hidden);
void set_highlighted(gui::ogl_object* object, bool highlighted);

// light properties
void add_default_lights(gui::ogl_scene* scene);
void set_light(gui::ogl_light* light, const vec3f& position,
    const vec3f& emission, ogl_light_type type, bool camera);

// light size
void clear_lights(gui::ogl_scene* scene);
bool has_max_lights(gui::ogl_scene* scene);

// Draw an OpenGL scene
void draw_scene(gui::ogl_scene* scene, gui::ogl_camera* camera,
    const vec4i& viewport, const ogl_scene_params& params);

}  // namespace yocto::gui

// -----------------------------------------------------------------------------
// UI APPLICATION
// -----------------------------------------------------------------------------
namespace yocto::gui {

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
using init_callback =
    std::function<void(gui::gui_window*, const gui_input& input)>;
// Clear callback called after the window is cloased
using clear_callback =
    std::function<void(gui::gui_window*, const gui_input& input)>;
// Draw callback called every frame and when resizing
using draw_callback =
    std::function<void(gui::gui_window*, const gui_input& input)>;
// Draw callback for drawing widgets
using widgets_callback =
    std::function<void(gui::gui_window*, const gui_input& input)>;
// Drop callback that returns that list of dropped strings.
using drop_callback = std::function<void(
    gui_window*, const std::vector<std::string>&, const gui_input& input)>;
// Key callback that returns key codes, pressed/released flag and modifier keys
using key_callback = std::function<void(
    gui::gui_window*, int key, bool pressed, const gui_input& input)>;
// Char callback that returns ASCII key
using char_callback = std::function<void(
    gui::gui_window*, unsigned int key, const gui_input& input)>;
// Mouse click callback that returns left/right button, pressed/released flag,
// modifier keys
using click_callback = std::function<void(
    gui::gui_window*, bool left, bool pressed, const gui_input& input)>;
// Scroll callback that returns scroll amount
using scroll_callback =
    std::function<void(gui::gui_window*, float amount, const gui_input& input)>;
// Update functions called every frame
using uiupdate_callback =
    std::function<void(gui::gui_window*, const gui_input& input)>;
// Update functions called every frame
using update_callback =
    std::function<void(gui::gui_window*, const gui_input& input)>;

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
void run_ui(const vec2i& size, const std::string& title,
    const gui_callbacks& callbaks, int widgets_width = 320,
    bool widgets_left = true);

}  // namespace yocto::gui

// -----------------------------------------------------------------------------
// UI WINDOW
// -----------------------------------------------------------------------------
namespace yocto::gui {

// OpenGL window wrapper
struct gui_window {
  GLFWwindow*   win           = nullptr;
  std::string       title         = "";
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
void init_window(gui::gui_window* win, const vec2i& size,
    const std::string& title, bool widgets, int widgets_width = 320,
    bool widgets_left = true);

// Window cleanup
void clear_window(gui::gui_window* win);

// Set callbacks
void set_init_callback(gui::gui_window* win, init_callback init_cb);
void set_clear_callback(gui::gui_window* win, clear_callback clear_cb);
void set_draw_callback(gui::gui_window* win, draw_callback draw_cb);
void set_widgets_callback(gui::gui_window* win, widgets_callback widgets_cb);
void set_drop_callback(gui::gui_window* win, drop_callback drop_cb);
void set_key_callback(gui::gui_window* win, key_callback cb);
void set_char_callback(gui::gui_window* win, char_callback cb);
void set_click_callback(gui::gui_window* win, click_callback cb);
void set_scroll_callback(gui::gui_window* win, scroll_callback cb);
void set_uiupdate_callback(gui::gui_window* win, uiupdate_callback cb);
void set_update_callback(gui::gui_window* win, update_callback cb);

// Run loop
void run_ui(gui::gui_window* win);
void set_close(gui::gui_window* win, bool close);

}  // namespace yocto::gui

// -----------------------------------------------------------------------------
// OPENGL WIDGETS
// -----------------------------------------------------------------------------
namespace yocto::gui {

bool begin_header(gui::gui_window* win, const char* title);
void end_header(gui::gui_window* win);

void draw_label(gui::gui_window* win, const char* lbl, const std::string& text);

void draw_separator(gui::gui_window* win);
void continue_line(gui::gui_window* win);

bool draw_button(gui::gui_window* win, const char* lbl, bool enabled = true);

bool draw_textinput(gui::gui_window* win, const char* lbl, std::string& value);

bool draw_slider(
    gui_window* win, const char* lbl, float& value, float min, float max);
bool draw_slider(
    gui_window* win, const char* lbl, vec2f& value, float min, float max);
bool draw_slider(
    gui_window* win, const char* lbl, vec3f& value, float min, float max);
bool draw_slider(
    gui_window* win, const char* lbl, vec4f& value, float min, float max);

bool draw_slider(
    gui::gui_window* win, const char* lbl, int& value, int min, int max);
bool draw_slider(
    gui::gui_window* win, const char* lbl, vec2i& value, int min, int max);
bool draw_slider(
    gui::gui_window* win, const char* lbl, vec3i& value, int min, int max);
bool draw_slider(
    gui::gui_window* win, const char* lbl, vec4i& value, int min, int max);

bool draw_dragger(gui::gui_window* win, const char* lbl, float& value,
    float speed = 1.0f, float min = 0.0f, float max = 0.0f);
bool draw_dragger(gui::gui_window* win, const char* lbl, vec2f& value,
    float speed = 1.0f, float min = 0.0f, float max = 0.0f);
bool draw_dragger(gui::gui_window* win, const char* lbl, vec3f& value,
    float speed = 1.0f, float min = 0.0f, float max = 0.0f);
bool draw_dragger(gui::gui_window* win, const char* lbl, vec4f& value,
    float speed = 1.0f, float min = 0.0f, float max = 0.0f);

bool draw_dragger(gui::gui_window* win, const char* lbl, int& value,
    float speed = 1, int min = 0, int max = 0);
bool draw_dragger(gui::gui_window* win, const char* lbl, vec2i& value,
    float speed = 1, int min = 0, int max = 0);
bool draw_dragger(gui::gui_window* win, const char* lbl, vec3i& value,
    float speed = 1, int min = 0, int max = 0);
bool draw_dragger(gui::gui_window* win, const char* lbl, vec4i& value,
    float speed = 1, int min = 0, int max = 0);

bool draw_checkbox(gui::gui_window* win, const char* lbl, bool& value);

bool draw_coloredit(gui::gui_window* win, const char* lbl, vec3f& value);
bool draw_coloredit(gui::gui_window* win, const char* lbl, vec4f& value);

bool draw_hdrcoloredit(gui::gui_window* win, const char* lbl, vec3f& value);
bool draw_hdrcoloredit(gui::gui_window* win, const char* lbl, vec4f& value);

bool draw_combobox(gui::gui_window* win, const char* lbl, int& idx,
    const std::vector<std::string>& labels);
bool draw_combobox(gui::gui_window* win, const char* lbl, std::string& value,
    const std::vector<std::string>& labels);
bool draw_combobox(gui::gui_window* win, const char* lbl, int& idx, int num,
    const std::function<std::string(int)>& labels, bool include_null = false);

template <typename T>
inline bool draw_combobox(gui::gui_window* win, const char* lbl, T*& value,
    const std::vector<T*>& vals, bool include_null = false) {
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
inline bool draw_combobox(gui::gui_window* win, const char* lbl, T*& value,
    const std::vector<T*>& vals, const std::vector<std::string>& labels,
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

void draw_progressbar(gui::gui_window* win, const char* lbl, float fraction);
void draw_progressbar(
    gui::gui_window* win, const char* lbl, int current, int total);

void draw_histogram(
    gui::gui_window* win, const char* lbl, const std::vector<float>& values);
void draw_histogram(
    gui::gui_window* win, const char* lbl, const std::vector<vec2f>& values);
void draw_histogram(
    gui::gui_window* win, const char* lbl, const std::vector<vec3f>& values);
void draw_histogram(
    gui::gui_window* win, const char* lbl, const std::vector<vec4f>& values);

bool draw_messages(gui::gui_window* win);
void push_message(gui::gui_window* win, const std::string& message);
bool draw_filedialog(gui::gui_window* win, const char* lbl, std::string& path,
    bool save, const std::string& dirname, const std::string& filename,
    const std::string& filter);
bool draw_filedialog_button(gui::gui_window* win, const char* button_lbl,
    bool button_active, const char* lbl, std::string& path, bool save,
    const std::string& dirname, const std::string& filename,
    const std::string& filter);

void log_info(gui::gui_window* win, const std::string& msg);
void log_error(gui::gui_window* win, const std::string& msg);
void clear_log(gui::gui_window* win);
void draw_log(gui::gui_window* win);

}  // namespace yocto::gui

#endif
