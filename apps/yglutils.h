//
// Utilities to use OpenGL 3, GLFW and ImGui.
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

#ifndef _YGLUTILS_H_
#define _YGLUTILS_H_

#include "../yocto/ygl.h"

// forward declaration
struct GLFWwindow;

namespace ygl {

void clear_glframebuffer(const vec4f& color, bool clear_depth = true);

void set_glviewport(int x, int y, int w, int h);
void set_glviewport(const vec2i& size);

void set_glwireframe(bool enabled);
void set_glblending(bool enabled);

struct glprogram {
    uint program_id = 0;
    uint vertex_shader_id = 0;
    uint fragment_shader_id = 0;
    uint vertex_array_object_id = 0;

         operator bool() const { return (bool)program_id; }
};

bool init_glprogram(
    glprogram& program, const char* vertex, const char* fragment);

void delete_glprogram(glprogram& program);

void bind_glprogram(glprogram& program);
void unbind_glprogram();

struct gltexture {
    uint texture_id    = 0;
    int  width  = 0;
    int  height = 0;

    operator bool() const { return (bool)texture_id; }
};

bool init_gltexture(gltexture& texture, const image<vec4f>& img, bool as_float,
    bool linear, bool mipmap);
void update_gltexture(gltexture& texture, const image<vec4f>& img,
    bool as_float, bool linear, bool mipmap);

bool init_gltexture(gltexture& texture, const image<vec4b>& img, bool as_srgb,
    bool linear, bool mipmap);
void update_gltexture(gltexture& texture, const image<vec4b>& img, bool as_srgb,
    bool linear, bool mipmap);

void delete_gltexture(gltexture& texture);

struct glarraybuffer {
    uint buffer_id       = 0;
    int  num       = 0;
    int  elem_size = 0;

    operator bool() const { return (bool)buffer_id; }
};

struct glelementbuffer {
    uint buffer_id       = 0;
    int  num       = 0;
    int  elem_size = 0;

    operator bool() const { return (bool)buffer_id; }
};

bool init_glarraybuffer(
    glarraybuffer& buffer, const vector<float>& data, bool dynamic = false);
bool init_glarraybuffer(
    glarraybuffer& buffer, const vector<vec2f>& data, bool dynamic = false);
bool init_glarraybuffer(
    glarraybuffer& buffer, const vector<vec3f>& data, bool dynamic = false);
bool init_glarraybuffer(
    glarraybuffer& buffer, const vector<vec4f>& data, bool dynamic = false);

void delete_glarraybuffer(glarraybuffer& buffer);

bool init_glelementbuffer(
    glelementbuffer& buffer, const vector<int>& data, bool dynamic = false);
bool init_glelementbuffer(
    glelementbuffer& buffer, const vector<vec2i>& data, bool dynamic = false);
bool init_glelementbuffer(
    glelementbuffer& buffer, const vector<vec3i>& data, bool dynamic = false);

void delete_glelementbuffer(glelementbuffer& buffer);

int get_gluniform_location(const glprogram& program, const char* name);

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
    const glprogram& program, const char* name, const T& value) {
    set_gluniform(get_gluniform_location(program, name), value);
}

void set_gluniform_texture(int locatiom, const gltexture& texture, int unit);
void set_gluniform_texture(
    glprogram& program, const char* name, const gltexture& texture, int unit);
void set_gluniform_texture(
    int locatiom, int locatiom_on, const gltexture& texture, int unit);
void set_gluniform_texture(glprogram& program, const char* name,
    const char* name_on, const gltexture& texture, int unit);

int get_glvertexattrib_location(const glprogram& program, const char* name);

void set_glvertexattrib(int locatiom, const glarraybuffer& buffer, float value);
void set_glvertexattrib(
    int locatiom, const glarraybuffer& buffer, const vec2f& value);
void set_glvertexattrib(
    int locatiom, const glarraybuffer& buffer, const vec3f& value);
void set_glvertexattrib(
    int locatiom, const glarraybuffer& buffer, const vec4f& value);

template <typename T>
inline void set_glvertexattrib(const glprogram& program, const char* name,
    const glarraybuffer& buffer, const T& value) {
    set_glvertexattrib(get_glvertexattrib_location(program, name), buffer, value);
}

void draw_glpoints(const glelementbuffer& buffer, int num);
void draw_gllines(const glelementbuffer& buffer, int num);
void draw_gltriangles(const glelementbuffer& buffer, int num);

void draw_glimage(const gltexture& texture, const vec2i& image_size,
    const vec2i& window_size, const vec2f& image_center, float image_scale);
void draw_glimage_background(const vec2i& image_size, const vec2i& window_size,
    const vec2f& image_center, float image_scale, float border_size = 2);

struct glwindow;
using refresh_glcallback = function<void(const glwindow&)>;
using drop_glcallback = function<void(const glwindow&, const vector<string>&)>;

struct glwindow {
    GLFWwindow*        win        = nullptr;
    void*              user_ptr   = nullptr;
    refresh_glcallback refresh_cb = {};
    drop_glcallback    drop_cb    = {};
};

bool init_glwindow(glwindow& win, int width, int height, const string& title,
    void* user_pointer, refresh_glcallback refresh_cb);
void delete_glwindow(glwindow& win);

void set_drop_glcallback(glwindow& win, drop_glcallback drop_cb);

void* get_user_pointer(const glwindow& win);

vec2i get_glframebuffer_size(const glwindow& win);
vec2i get_glwindow_size(const glwindow& win);

bool should_glwindow_close(const glwindow& win);

vec2f get_glmouse_pos(const glwindow& win);
bool  get_glmouse_left(const glwindow& win);
bool  get_glmouse_right(const glwindow& win);
bool  get_glalt_key(const glwindow& win);
bool  get_glshift_key(const glwindow& win);

void process_glevents(const glwindow& win, bool wait = false);
void swap_glbuffers(const glwindow& win);

void init_glwidgets(const glwindow& win);
bool get_glwidgets_active(const glwindow& win);

void begin_glwidgets_frame(const glwindow& win);
void end_glwidgets_frame(const glwindow& win);

bool begin_glwidgets_window(const glwindow& win, const char* title);

bool begin_header_glwidget(const glwindow& win, const char* title);
void end_header_glwidget(const glwindow& win);

void draw_label_glwidgets(
    const glwindow& win, const char* lbl, const string& texture);
void draw_label_glwidgets(
    const glwindow& win, const char* lbl, const char* fmt, ...);

bool begin_header_widget(const glwindow& win, const char* label);
void end_header_widget(const glwindow& win);

void draw_separator_glwidget(const glwindow& win);
void continue_glwidgets_line(const glwindow& win);

bool draw_button_glwidget(const glwindow& win, const char* lbl);

bool draw_textinput_glwidget(const glwindow& win, const char* lbl, string& value);
bool draw_slider_glwidget(
    const glwindow& win, const char* lbl, float& value, float min, float max);
bool draw_slider_glwidget(
    const glwindow& win, const char* lbl, vec2f& value, float min, float max);
bool draw_slider_glwidget(
    const glwindow& win, const char* lbl, vec3f& value, float min, float max);
bool draw_slider_glwidget(
    const glwindow& win, const char* lbl, vec4f& value, float min, float max);

bool draw_slider_glwidget(
    const glwindow& win, const char* lbl, int& value, int min, int max);
bool draw_slider_glwidget(
    const glwindow& win, const char* lbl, vec2i& value, int min, int max);
bool draw_slider_glwidget(
    const glwindow& win, const char* lbl, vec3i& value, int min, int max);
bool draw_slider_glwidget(
    const glwindow& win, const char* lbl, vec4i& value, int min, int max);

bool draw_dragger_glwidget(const glwindow& win, const char* lbl, float& value,
    float speed = 1.0f, float min = 0.0f, float max = 0.0f);
bool draw_dragger_glwidget(const glwindow& win, const char* lbl, vec2f& value,
    float speed = 1.0f, float min = 0.0f, float max = 0.0f);
bool draw_dragger_glwidget(const glwindow& win, const char* lbl, vec3f& value,
    float speed = 1.0f, float min = 0.0f, float max = 0.0f);
bool draw_dragger_glwidget(const glwindow& win, const char* lbl, vec4f& value,
    float speed = 1.0f, float min = 0.0f, float max = 0.0f);

bool draw_dragger_glwidget(const glwindow& win, const char* lbl, int& value,
    float speed = 1, int min = 0, int max = 0);
bool draw_dragger_glwidget(const glwindow& win, const char* lbl, vec2i& value,
    float speed = 1, int min = 0, int max = 0);
bool draw_dragger_glwidget(const glwindow& win, const char* lbl, vec3i& value,
    float speed = 1, int min = 0, int max = 0);
bool draw_dragger_glwidget(const glwindow& win, const char* lbl, vec4i& value,
    float speed = 1, int min = 0, int max = 0);

bool draw_checkbox_glwidget(const glwindow& win, const char* lbl, bool& value);

bool draw_coloredit_glwidget(const glwindow& win, const char* lbl, vec3f& value);
bool draw_coloredit_glwidget(const glwindow& win, const char* lbl, vec4f& value);

bool begin_treenode_glwidget(const glwindow& win, const char* lbl);
void end_treenode_glwidget(const glwindow& win);

bool begin_selectabletreenode_glwidget(
    const glwindow& win, const char* lbl, bool& selected);
void begin_selectabletreeleaf_glwidget(
    const glwindow& win, const char* lbl, bool& selected);

bool draw_combobox_glwidget(const glwindow& win, const char* lbl, int& idx,
    const vector<string>& labels);
bool draw_combobox_glwidget(const glwindow& win, const char* lbl, string& value,
    const vector<string>& labels);
bool draw_combobox_glwidget(const glwindow& win, const char* lbl, int& idx,
    int num, const function<const char*(int)>& labels, bool include_null = false);

template <typename T>
inline bool draw_combobox_glwidget(const glwindow& win, const char* lbl,
    int& idx, const vector<T*>& vals, bool include_null = false) {
    return draw_combobox_glwidget(win, lbl, idx, (int)vals.size(),
        [&](int idx) { return vals[idx]->name.c_str(); }, include_null);
}
template <typename T>
inline bool draw_combobox_glwidget(const glwindow& win, const char* lbl,
    int& idx, const vector<T>& vals, bool include_null = false) {
    return draw_combobox_glwidget(win, lbl, idx, (int)vals.size(),
        [&](int idx) { return vals[idx].name.c_str(); }, include_null);
}
template <typename T>
inline bool draw_combobox_glwidget(const glwindow& win, const char* lbl,
    int& idx, const deque<T>& vals, bool include_null = false) {
    return draw_combobox_glwidget(win, lbl, idx, (int)vals.size(),
        [&](int idx) { return vals[idx].name.c_str(); }, include_null);
}

void begin_child_glwidget(
    const glwindow& win, const char* lbl, const vec2i& size);
void end_child_glwidget(const glwindow& win);

bool begin_popup_modal(const glwindow& win, const char* lbl);
void end_popup_modal(const glwindow& win);

}  // namespace ygl

#endif
