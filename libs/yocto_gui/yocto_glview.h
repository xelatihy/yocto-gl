//
// Yocto/GlView: Simple OpenGL viewers.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
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

#ifndef _YOCTO_GLVIEW_
#define _YOCTO_GLVIEW_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <yocto/yocto_image.h>
#include <yocto/yocto_scene.h>
#include <yocto/yocto_trace.h>

#include <array>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::string;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE AND TRACE VIEW
// -----------------------------------------------------------------------------
namespace yocto {

// Open a window and show an image
void view_image(
    const string& title, const string& name, const color_image& image);

// Open a window and show a set of images
void view_images(const string& title, const vector<string>& names,
    const vector<color_image>& images);

// Open a window and show an image for color grading
void colorgrade_image(
    const string& title, const string& name, const color_image& image);

// Open a window and show an scene via path tracing
void view_scene(const string& title, const string& name, scene_model& scene,
    const trace_params& params = {}, bool print = true, bool edit = false);

// GUI callback
struct glinput_state;
using glview_callback = std::function<void(const glinput_state& input,
    vector<int>& updated_shapes, vector<int>& updated_textures)>;

// Open a window and show an scene via OpenGL shading
struct glscene_params;
void glview_scene(const string& title, const string& name, scene_model& scene,
    const glscene_params& params, const glview_callback& widgets_callback = {},
    const glview_callback& uiupdate_callback = {},
    const glview_callback& update_callback   = {});

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE DRAWING
// -----------------------------------------------------------------------------
namespace yocto {

// OpenGL image data
struct glimage_state {
  // image properties
  int width  = 0;
  int height = 0;

  // Opengl state
  uint texture     = 0;  // texture
  uint program     = 0;  // program
  uint vertex      = 0;
  uint fragment    = 0;
  uint vertexarray = 0;  // vertex
  uint positions   = 0;
  uint triangles   = 0;  // elements
};

// create image drawing program
bool init_image(glimage_state& glimage);

// clear image
void clear_image(glimage_state& glimage);

// update image data
void set_image(glimage_state& glimage, const color_image& image);

// OpenGL image drawing params
struct glimage_params {
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
void draw_image(glimage_state& image, const glimage_params& params);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE DRAWING
// -----------------------------------------------------------------------------
namespace yocto {

// Opengl texture
struct glscene_texture {
  // texture properties
  int width  = 0;
  int height = 0;

  // opengl state
  uint texture = 0;
};

// Create texture
void set_texture(glscene_texture& gltexture, const scene_texture& texture);

// Clean texture
void clear_texture(glscene_texture& gltexture);

// Opengl shape
struct glscene_shape {
  // Shape properties
  int num_positions = 0;
  int num_normals   = 0;
  int num_texcoords = 0;
  int num_colors    = 0;
  int num_tangents  = 0;
  int num_points    = 0;
  int num_lines     = 0;
  int num_triangles = 0;
  int num_quads     = 0;

  // OpenGl state
  uint  vertexarray = 0;
  uint  positions   = 0;
  uint  normals     = 0;
  uint  texcoords   = 0;
  uint  colors      = 0;
  uint  tangents    = 0;
  uint  points      = 0;
  uint  lines       = 0;
  uint  triangles   = 0;
  uint  quads       = 0;
  float point_size  = 1;
};

// Create shape
void set_shape(glscene_shape& glshape, const scene_shape& shape);

// Clean shape
void clear_shape(glscene_shape& glshape);

// Opengl scene
struct glscene_state {
  // scene objects
  vector<glscene_shape>   shapes   = {};
  vector<glscene_texture> textures = {};

  // programs
  uint program  = 0;
  uint vertex   = 0;
  uint fragment = 0;
};

// Shading type
enum struct glscene_lighting_type { camlight, eyelight };

// Shading name
const auto glscene_lighting_names = vector<string>{"camlight", "eyelight"};

// Draw options
struct glscene_params {
  int                   camera           = 0;
  int                   resolution       = 1280;
  bool                  wireframe        = false;
  glscene_lighting_type lighting         = glscene_lighting_type::camlight;
  float                 exposure         = 0;
  float                 gamma            = 2.2f;
  bool                  faceted          = false;
  bool                  double_sided     = true;
  bool                  non_rigid_frames = true;
  float                 near             = 0.01f;
  float                 far              = 10000.0f;
  bool                  hide_environment = false;
  vec4f                 background       = vec4f{0.15f, 0.15f, 0.15f, 1.0f};
};

// init scene
void init_glscene(glscene_state& glscene, const scene_model& scene);

// update scene
void update_glscene(glscene_state& glscene, const scene_model& scene,
    const vector<int>& updated_shapes, const vector<int>& updated_textures);

// Clear an OpenGL scene
void clear_scene(glscene_state& scene);

// draw scene
void draw_scene(glscene_state& glscene, const scene_model& scene,
    const vec4i& viewport, const glscene_params& params);

}  // namespace yocto

// -----------------------------------------------------------------------------
// WINDOW
// -----------------------------------------------------------------------------
namespace yocto {

// Input state
struct glinput_state {
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
using init_glcallback = function<void(const glinput_state& input)>;
// Clear callback called after the window is cloased
using clear_glcallback = function<void(const glinput_state& input)>;
// Draw callback called every frame and when resizing
using draw_glcallback = function<void(const glinput_state& input)>;
// Draw callback for drawing widgets
using widgets_glcallback = function<void(const glinput_state& input)>;
// Update functions called every frame
using update_glcallback = function<void(const glinput_state& input)>;
// Update functions called every frame
using uiupdate_glcallback = function<void(const glinput_state& input)>;

// User interface callcaks
struct glwindow_callbacks {
  init_glcallback     init_cb     = {};
  clear_glcallback    clear_cb    = {};
  draw_glcallback     draw_cb     = {};
  widgets_glcallback  widgets_cb  = {};
  update_glcallback   update_cb   = {};
  uiupdate_glcallback uiupdate_cb = {};
};

// run the user interface with the give callbacks
void run_ui(const vec2i& size, const string& title,
    const glwindow_callbacks& callbaks, int widgets_width = 320,
    bool widgets_left = true);

}  // namespace yocto

// -----------------------------------------------------------------------------
// WIDGETS
// -----------------------------------------------------------------------------
namespace yocto {

// Headers
bool begin_glheader(const char* title);
void end_glheader();

// Labels
void draw_gllabel(const char* lbl, const string& text);
void draw_gllabel(const char* lbl, int value);
void draw_gllabel(const char* lbl, bool value);

// Lines
void draw_glseparator();
void continue_glline();

// Buttons
bool draw_glbutton(const char* lbl, bool enabled = true);

// Text
bool draw_gltextinput(const char* lbl, string& value);

// Slider
bool draw_glslider(const char* lbl, float& value, float min, float max);
bool draw_glslider(const char* lbl, vec2f& value, float min, float max);
bool draw_glslider(const char* lbl, vec3f& value, float min, float max);
bool draw_glslider(const char* lbl, vec4f& value, float min, float max);
bool draw_glslider(const char* lbl, int& value, int min, int max);
bool draw_glslider(const char* lbl, vec2i& value, int min, int max);
bool draw_glslider(const char* lbl, vec3i& value, int min, int max);
bool draw_glslider(const char* lbl, vec4i& value, int min, int max);

// Dragger
bool draw_gldragger(const char* lbl, float& value, float speed = 1.0f,
    float min = 0.0f, float max = 0.0f);
bool draw_gldragger(const char* lbl, vec2f& value, float speed = 1.0f,
    float min = 0.0f, float max = 0.0f);
bool draw_gldragger(const char* lbl, vec3f& value, float speed = 1.0f,
    float min = 0.0f, float max = 0.0f);
bool draw_gldragger(const char* lbl, vec4f& value, float speed = 1.0f,
    float min = 0.0f, float max = 0.0f);
bool draw_gldragger(
    const char* lbl, int& value, float speed = 1, int min = 0, int max = 0);
bool draw_gldragger(
    const char* lbl, vec2i& value, float speed = 1, int min = 0, int max = 0);
bool draw_gldragger(
    const char* lbl, vec3i& value, float speed = 1, int min = 0, int max = 0);
bool draw_gldragger(
    const char* lbl, vec4i& value, float speed = 1, int min = 0, int max = 0);

// Checkbox
bool draw_glcheckbox(const char* lbl, bool& value);
bool draw_glcheckbox(const char* lbl, bool& value, bool invert);

// Color editor
bool draw_glcoloredit(const char* lbl, vec3f& value);
bool draw_glcoloredit(const char* lbl, vec4f& value);
bool draw_glcoloredit(const char* lbl, vec4b& value);
bool draw_glcoloredithdr(const char* lbl, vec3f& value);
bool draw_glcoloredithdr(const char* lbl, vec4f& value);

// Combo box
bool draw_glcombobox(const char* lbl, int& idx, const vector<string>& labels,
    bool include_null = false);
bool draw_glcombobox(const char* lbl, string& value,
    const vector<string>& labels, bool include_null = false);

// Progress bar
void draw_glprogressbar(const char* lbl, float fraction);
void draw_glprogressbar(const char* lbl, int current, int total);

}  // namespace yocto

#endif
