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

#include "../yocto/yocto_scene.h"

#include <functional>

// forward declaration
struct GLFWwindow;

// -----------------------------------------------------------------------------
// OPENGL FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

void clear_glframebuffer(const vec4f& color, bool clear_depth = true);

void set_glviewport(const vec4i& viewport);

void set_glwireframe(bool enabled);
void set_glblending(bool enabled);

struct opengl_program {
  uint program_id             = 0;
  uint vertex_shader_id       = 0;
  uint fragment_shader_id     = 0;
  uint vertex_array_object_id = 0;

  operator bool() const { return (bool)program_id; }
};

void init_glprogram(
    opengl_program& program, const char* vertex, const char* fragment);

void delete_glprogram(opengl_program& program);

void bind_glprogram(opengl_program& program);
void unbind_opengl_program();

struct opengl_texture {
  uint  texture_id = 0;
  vec2i size       = {0, 0};

  operator bool() const { return (bool)texture_id; }
};

void init_gltexture(opengl_texture& texture, const vec2i& size, bool as_float,
    bool as_srgb, bool linear, bool mipmap);

void update_gltexture(
    opengl_texture& texture, const image<vec4f>& img, bool mipmap);
void update_gltexture_region(opengl_texture& texture, const image<vec4f>& img,
    const image_region& region, bool mipmap);

inline void init_gltexture(opengl_texture& texture, const image<vec4f>& img,
    bool as_float, bool linear, bool mipmap) {
  init_gltexture(texture, img.size(), as_float, false, linear, mipmap);
  update_gltexture(texture, img, mipmap);
}

void init_gltexture(opengl_texture& texture, const image<vec4b>& img,
    bool as_srgb, bool linear, bool mipmap);
void update_gltexture(
    opengl_texture& texture, const image<vec4b>& img, bool mipmap);
void update_gltexture_region(opengl_texture& texture, const image<vec4b>& img,
    const image_region& region, bool mipmap);

inline void init_gltexture(opengl_texture& texture, const image<vec4b>& img,
    bool as_srgb, bool linear, bool mipmap) {
  init_gltexture(texture, img.size(), false, as_srgb, linear, mipmap);
  update_gltexture(texture, img, mipmap);
}

void delete_gltexture(opengl_texture& texture);

struct opengl_arraybuffer {
  uint buffer_id = 0;
  int  num       = 0;
  int  elem_size = 0;

  operator bool() const { return (bool)buffer_id; }
};

struct opengl_elementbuffer {
  uint buffer_id = 0;
  int  num       = 0;
  int  elem_size = 0;

  operator bool() const { return (bool)buffer_id; }
};

void init_glarraybuffer(opengl_arraybuffer& buffer, const vector<float>& data,
    bool dynamic = false);
void init_glarraybuffer(opengl_arraybuffer& buffer, const vector<vec2f>& data,
    bool dynamic = false);
void init_glarraybuffer(opengl_arraybuffer& buffer, const vector<vec3f>& data,
    bool dynamic = false);
void init_glarraybuffer(opengl_arraybuffer& buffer, const vector<vec4f>& data,
    bool dynamic = false);

void delete_glarraybuffer(opengl_arraybuffer& buffer);

void init_glelementbuffer(opengl_elementbuffer& buffer, const vector<int>& data,
    bool dynamic = false);
void init_glelementbuffer(opengl_elementbuffer& buffer,
    const vector<vec2i>& data, bool dynamic = false);
void init_glelementbuffer(opengl_elementbuffer& buffer,
    const vector<vec3i>& data, bool dynamic = false);

void delete_glelementbuffer(opengl_elementbuffer& buffer);

int get_gluniform_location(const opengl_program& program, const char* name);

void set_gluniform(int locatiom, int value);
void set_gluniform(int locatiom, const vec2i& value);
void set_gluniform(int locatiom, const vec3i& value);
void set_gluniform(int locatiom, const vec4i& value);
void set_gluniform(int locatiom, float value);
void set_gluniform(int locatiom, const vec2f& value);
void set_gluniform(int locatiom, const vec3f& value);
void set_gluniform(int locatiom, const vec4f& value);
void set_gluniform(int locatiom, const mat4f& value);
void set_gluniform(int locatiom, const frame3f& value);

template <typename T>
inline void set_gluniform(
    const opengl_program& program, const char* name, const T& value) {
  set_gluniform(get_gluniform_location(program, name), value);
}

void set_gluniform_texture(
    int locatiom, const opengl_texture& texture, int unit);
void set_gluniform_texture(opengl_program& program, const char* name,
    const opengl_texture& texture, int unit);
void set_gluniform_texture(
    int locatiom, int locatiom_on, const opengl_texture& texture, int unit);
void set_gluniform_texture(opengl_program& program, const char* name,
    const char* name_on, const opengl_texture& texture, int unit);

int get_glvertexattrib_location(
    const opengl_program& program, const char* name);

void set_glvertexattrib(
    int locatiom, const opengl_arraybuffer& buffer, float value);
void set_glvertexattrib(
    int locatiom, const opengl_arraybuffer& buffer, const vec2f& value);
void set_glvertexattrib(
    int locatiom, const opengl_arraybuffer& buffer, const vec3f& value);
void set_glvertexattrib(
    int locatiom, const opengl_arraybuffer& buffer, const vec4f& value);

template <typename T>
inline void set_glvertexattrib(const opengl_program& program, const char* name,
    const opengl_arraybuffer& buffer, const T& value) {
  set_glvertexattrib(get_glvertexattrib_location(program, name), buffer, value);
}

void draw_glpoints(const opengl_elementbuffer& buffer, int num);
void draw_gllines(const opengl_elementbuffer& buffer, int num);
void draw_gltriangles(const opengl_elementbuffer& buffer, int num);

void draw_glimage(const opengl_texture& texture, int win_width, int win_height,
    const vec2f& image_center, float image_scale);
void draw_glimage_background(const opengl_texture& texture, int win_width,
    int win_height, const vec2f& image_center, float image_scale,
    float border_size = 2);

struct opengl_window;
using refresh_glcallback = std::function<void(const opengl_window&)>;
using drop_glcallback =
    std::function<void(const opengl_window&, const vector<string>&)>;

struct opengl_window {
  GLFWwindow*        win           = nullptr;
  void*              user_ptr      = nullptr;
  refresh_glcallback refresh_cb    = {};
  drop_glcallback    drop_cb       = {};
  int                widgets_width = 0;
  bool               widgets_left  = true;
};

void init_glwindow(opengl_window& win, const vec2i& size, const string& title,
    void* user_pointer, refresh_glcallback refresh_cb);
void delete_glwindow(opengl_window& win);

void set_drop_glcallback(opengl_window& win, drop_glcallback drop_cb);

void* get_gluser_pointer(const opengl_window& win);

vec2i get_glwindow_size(const opengl_window& win, bool ignore_widgets = true);
vec2i get_glframebuffer_size(
    const opengl_window& win, bool ignore_widgets = true);
vec4i get_glframebuffer_viewport(
    const opengl_window& win, bool ignore_widgets = true);

bool should_glwindow_close(const opengl_window& win);
void set_glwindow_close(const opengl_window& win, bool close);

vec2f get_glmouse_pos(const opengl_window& win, bool ignore_widgets = true);
bool  get_glmouse_left(const opengl_window& win);
bool  get_glmouse_right(const opengl_window& win);
bool  get_glalt_key(const opengl_window& win);
bool  get_glshift_key(const opengl_window& win);

void process_glevents(const opengl_window& win, bool wait = false);
void swap_glbuffers(const opengl_window& win);

void init_glwidgets(opengl_window& win, int width = 320, bool left = true);
bool get_glwidgets_active(const opengl_window& win);

void begin_glwidgets(const opengl_window& win);
void end_glwidgets(const opengl_window& win);

bool begin_glwidgets_window(const opengl_window& win, const char* title);

bool begin_glheader(const opengl_window& win, const char* title);
void end_glheader(const opengl_window& win);

bool begin_gltabbar(const opengl_window& win, const char* title);
void end_gltabbar(const opengl_window& win);

bool begin_gltabitem(const opengl_window& win, const char* title);
void end_gltabitem(const opengl_window& win);

void open_glmodal(const opengl_window& win, const char* lbl);
void clear_glmodal(const opengl_window& win);
bool begin_glmodal(const opengl_window& win, const char* lbl);
void end_glmodal(const opengl_window& win);
bool is_glmodal_open(const opengl_window& win, const char* lbl);

bool draw_glmessages(const opengl_window& win);
void push_glmessage(const string& message);
void push_glmessage(const opengl_window& win, const string& message);
bool draw_glmessage(
    const opengl_window& win, const char* lbl, const string& message);
bool draw_glfiledialog(const opengl_window& win, const char* lbl, string& path,
    bool save, const string& dirname, const string& filename,
    const string& filter);
bool draw_glfiledialog_button(const opengl_window& win, const char* button_lbl,
    bool button_active, const char* lbl,
    string& path, bool save, const string& dirname, const string& filename,
    const string& filter);

void draw_gltext(const opengl_window& win, const string& text);
void draw_gllabel(
    const opengl_window& win, const char* lbl, const string& text);
void draw_gllabel(
    const opengl_window& win, const char* lbl, const char* fmt, ...);

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

bool begin_gltreenode(const opengl_window& win, const char* lbl);
void end_gltreenode(const opengl_window& win);

bool begin_glselectabletreenode(
    const opengl_window& win, const char* lbl, bool& selected);
void begin_glselectabletreeleaf(
    const opengl_window& win, const char* lbl, bool& selected);

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

void begin_glchild(
    const opengl_window& win, const char* lbl, const vec2i& size);
void end_glchild(const opengl_window& win);

void draw_glhistogram(
    const opengl_window& win, const char* lbl, const float* values, int count);
void draw_glhistogram(
    const opengl_window& win, const char* lbl, const vector<float>& values);
void draw_glhistogram(
    const opengl_window& win, const char* lbl, const vector<vec2f>& values);
void draw_glhistogram(
    const opengl_window& win, const char* lbl, const vector<vec3f>& values);
void draw_glhistogram(
    const opengl_window& win, const char* lbl, const vector<vec4f>& values);

void log_glinfo(const opengl_window& win, const char* msg);
void log_glinfo(const opengl_window& win, const string& msg);
void log_glerror(const opengl_window& win, const char* msg);
void log_glerror(const opengl_window& win, const string& msg);
void clear_gllogs(const opengl_window& win);
void draw_gllog(const opengl_window& win);

}  // namespace yocto

#endif
