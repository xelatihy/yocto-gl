//
// Yocto/GlView: Simple OpenGL viewers.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2022 Fabio Pellacini
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

#ifdef YOCTO_OPENGL

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <array>
#include <functional>
#include <string>
#include <vector>

#include "yocto_image.h"
#include "yocto_scene.h"
#include "yocto_trace.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::function;
using std::string;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE AND TRACE VIEW
// -----------------------------------------------------------------------------
namespace yocto {

// Open a window and show an image
void show_image_gui(
    const string& title, const string& name, const image_data& image);

// Open a window and show a set of images
void show_image_gui(const string& title, const vector<string>& names,
    const vector<image_data>& images);

// Open a window and show an image for color grading
void show_colorgrade_gui(
    const string& title, const string& name, const image_data& image);

// Open a window and show an scene via path tracing
void show_trace_gui(const string& title, const string& name, scene_data& scene,
    const trace_params& params = {}, bool print = true, bool edit = false);

// Open a window and show an scene via path tracing in cuda
void show_cutrace_gui(const string& title, const string& name,
    scene_data& scene, const trace_params& params = {}, bool print = true,
    bool edit = false);

// GUI callback
struct gui_input;
using glview_callback = std::function<void(const gui_input& input,
    vector<int>& updated_shapes, vector<int>& updated_textures)>;

// Shading type
enum struct shade_lighting { camlight, eyelight };

// Shading name
inline const auto shade_lighting_names = vector<string>{"camlight", "eyelight"};

// Draw options
struct shade_params {
  int            camera           = 0;
  int            resolution       = 1280;
  bool           wireframe        = false;
  shade_lighting lighting         = shade_lighting::camlight;
  float          exposure         = 0;
  float          gamma            = 2.2f;
  bool           faceted          = false;
  bool           double_sided     = true;
  bool           non_rigid_frames = true;
  float          near             = 0.01f;
  float          far              = 10000.0f;
  bool           hide_environment = false;
  vec4f          background       = vec4f{0.5f, 0.5f, 0.5f, 0.5f};
};

// Open a window and show an scene via OpenGL shading
struct shade_params;
void show_shade_gui(const string& title, const string& name, scene_data& scene,
    const shade_params& params, const glview_callback& widgets_callback = {},
    const glview_callback& uiupdate_callback = {},
    const glview_callback& update_callback   = {});

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE
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

// OpenGL image drawing params
struct glimage_params {
  vec2i window      = {512, 512};
  vec4i framebuffer = {0, 0, 512, 512};
  vec2f center      = {0, 0};
  float scale       = 1;
  bool  fit         = true;
  bool  checker     = true;
  float border_size = 2;
  vec4f background  = {0.5f, 0.5f, 0.5f, 1.0f};
  bool  tonemap     = false;
  float exposure    = 0;
  bool  srgb        = true;
  bool  filmic      = false;
};

// create image drawing program
bool init_image(glimage_state& glimage);

// clear image
void clear_image(glimage_state& glimage);

// update image data
void set_image(glimage_state& glimage, const image_data& image);

// draw image
void draw_image(glimage_state& image, const glimage_params& params);

}  // namespace yocto

// -----------------------------------------------------------------------------
// WINDOW
// -----------------------------------------------------------------------------
namespace yocto {

// Input state
struct gui_input {
  vec3i mouse       = {0, 0, 0};     // mouse buttons (left, right, middle)
  vec2f cursor      = {0, 0};        // position excluding widgets
  vec2f last        = {0, 0};        // last mouse position excluding widgets
  vec3i modifiers   = {0, 0, 0};     // modifieers (alt, shift, control)
  bool  onwidgets   = false;         // widgets are active
  vec2i window      = {0, 0};        // window size
  vec4i framebuffer = {0, 0, 0, 0};  // framebuffer viewport
};

// Init callback called after the window has opened
using gui_callback = function<void(const gui_input& input)>;

// User interface callcaks
struct gui_callbacks {
  gui_callback init     = {};
  gui_callback clear    = {};
  gui_callback draw     = {};
  gui_callback widgets  = {};
  gui_callback update   = {};
  gui_callback uiupdate = {};
};

// run the user interface with the give callbacks
void show_gui_window(const vec2i& size, const string& title,
    const gui_callbacks& callbaks, int widgets_width = 320,
    bool widgets_left = true);

}  // namespace yocto

// -----------------------------------------------------------------------------
// WIDGETS
// -----------------------------------------------------------------------------
namespace yocto {

// Headers
bool draw_gui_header(const char* title);
void end_gui_header();

// Labels
void draw_gui_label(const char* lbl, const string& text);
void draw_gui_label(const char* lbl, int value);
void draw_gui_label(const char* lbl, bool value);

// Lines
void draw_gui_separator();
void continue_gui_line();

// Buttons
bool draw_gui_button(const char* lbl, bool enabled = true);

// Text
bool draw_gui_textinput(const char* lbl, string& value);

// Slider
bool draw_gui_slider(const char* lbl, float& value, float min, float max);
bool draw_gui_slider(const char* lbl, vec2f& value, float min, float max);
bool draw_gui_slider(const char* lbl, vec3f& value, float min, float max);
bool draw_gui_slider(const char* lbl, vec4f& value, float min, float max);
bool draw_gui_slider(const char* lbl, int& value, int min, int max);
bool draw_gui_slider(const char* lbl, vec2i& value, int min, int max);
bool draw_gui_slider(const char* lbl, vec3i& value, int min, int max);
bool draw_gui_slider(const char* lbl, vec4i& value, int min, int max);

// Dragger
bool draw_gui_dragger(const char* lbl, float& value, float speed = 1.0f,
    float min = 0.0f, float max = 0.0f);
bool draw_gui_dragger(const char* lbl, vec2f& value, float speed = 1.0f,
    float min = 0.0f, float max = 0.0f);
bool draw_gui_dragger(const char* lbl, vec3f& value, float speed = 1.0f,
    float min = 0.0f, float max = 0.0f);
bool draw_gui_dragger(const char* lbl, vec4f& value, float speed = 1.0f,
    float min = 0.0f, float max = 0.0f);
bool draw_gui_dragger(
    const char* lbl, int& value, float speed = 1, int min = 0, int max = 0);
bool draw_gui_dragger(
    const char* lbl, vec2i& value, float speed = 1, int min = 0, int max = 0);
bool draw_gui_dragger(
    const char* lbl, vec3i& value, float speed = 1, int min = 0, int max = 0);
bool draw_gui_dragger(
    const char* lbl, vec4i& value, float speed = 1, int min = 0, int max = 0);

// Checkbox
bool draw_gui_checkbox(const char* lbl, bool& value);
bool draw_gui_checkbox(const char* lbl, bool& value, bool invert);

// Color editor
bool draw_gui_coloredit(const char* lbl, vec3f& value);
bool draw_gui_coloredit(const char* lbl, vec4f& value);
bool draw_gui_coloredit(const char* lbl, vec4b& value);
bool draw_gui_coloredithdr(const char* lbl, vec3f& value);
bool draw_gui_coloredithdr(const char* lbl, vec4f& value);

// Combo box
bool draw_gui_combobox(const char* lbl, int& idx, const vector<string>& labels,
    bool include_null = false);
bool draw_gui_combobox(const char* lbl, string& value,
    const vector<string>& labels, bool include_null = false);

// Progress bar
void draw_gui_progressbar(const char* lbl, float fraction);
void draw_gui_progressbar(const char* lbl, int current, int total);

}  // namespace yocto

// -----------------------------------------------------------------------------
// HIGH LEVEL WIDGETS
// -----------------------------------------------------------------------------
namespace yocto {

// draw tonemap params
bool draw_tonemap_widgets(
    const gui_input& input, float& exposure, bool& filmic);

// draw image inspector
bool draw_image_widgets(const gui_input& input, const image_data& image,
    const image_data& display, glimage_params& glparams);
bool draw_image_widgets(
    const gui_input& input, const image_data& image, glimage_params& glparams);

// update image params
void update_image_params(
    const gui_input& input, const image_data& image, glimage_params& glparams);

// update image params from mouse
bool uiupdate_image_params(const gui_input& input, glimage_params& glparams);

// update camera from mouse
bool uiupdate_camera_params(const gui_input& input, camera_data& camera);

// draw trace params
bool draw_trace_widgets(const gui_input& input, int sample,
    trace_params& params, const vector<string>& camera_names);

// scene selection
struct scene_selection {
  int camera      = 0;
  int instance    = 0;
  int environment = 0;
  int shape       = 0;
  int texture     = 0;
  int material    = 0;
  int subdiv      = 0;
};

// draw scene editor
bool draw_scene_widgets(scene_data& scene, scene_selection& selection,
    const function<void()>& before_edit = {});

}  // namespace yocto

#endif

#endif
