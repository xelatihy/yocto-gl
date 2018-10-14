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

struct glprogram {
    uint pid = 0;
    uint vid = 0;
    uint fid = 0;
    uint vao = 0;
         operator bool() const { return (bool)pid; }
};

glprogram make_glprogram(const char* vertex, const char* fragment);
void      bind_glprogram(glprogram& pid);
void      unbind_glprogram();

struct gltexture {
    uint tid = 0;
         operator bool() const { return (bool)tid; }
};

gltexture make_gltexture(
    const image<vec4f>& img, bool as_float, bool linear, bool mipmap);
void update_gltexture(gltexture& txt, const image<vec4f>& img, bool as_float,
    bool linear, bool mipmap);

gltexture make_gltexture(
    const image<vec4b>& img, bool as_srgb, bool linear, bool mipmap);
void update_gltexture(gltexture& txt, const image<vec4b>& img, bool as_srgb,
    bool linear, bool mipmap);

struct glarraybuffer {
    uint bid = 0;
         operator bool() const { return (bool)bid; }
};

struct glelementbuffer {
    uint bid = 0;
         operator bool() const { return (bool)bid; }
};

glarraybuffer make_glarraybuffer(const vector<float>& buf, bool dynamic = false);
glarraybuffer make_glarraybuffer(const vector<vec2f>& buf, bool dynamic = false);
glarraybuffer make_glarraybuffer(const vector<vec3f>& buf, bool dynamic = false);
glarraybuffer make_glarraybuffer(const vector<vec4f>& buf, bool dynamic = false);

glelementbuffer make_glelementbuffer(
    const vector<int>& buf, bool dynamic = false);
glelementbuffer make_glelementbuffer(
    const vector<vec2i>& buf, bool dynamic = false);
glelementbuffer make_glelementbuffer(
    const vector<vec3i>& buf, bool dynamic = false);

int get_gluniform_location(const glprogram& prog, const char* name);

void set_gluniform(int loc, int val);
void set_gluniform(int loc, const vec2i& val);
void set_gluniform(int loc, const vec3i& val);
void set_gluniform(int loc, const vec4i& val);
void set_gluniform(int loc, float val);
void set_gluniform(int loc, const vec2f& val);
void set_gluniform(int loc, const vec3f& val);
void set_gluniform(int loc, const vec4f& val);
void set_gluniform(int loc, const mat4f& val);
void set_gluniform(int loc, const frame3f& val);

template <typename T>
inline void set_gluniform(const glprogram& prog, const char* var, const T& val) {
    set_gluniform(get_gluniform_location(prog, var), val);
}

void set_gluniform_texture(int loc, const gltexture& txt, int unit);
void set_gluniform_texture(
    glprogram& prog, const char* var, const gltexture& txt, int unit);
void set_gluniform_texture(int loc, int loc_on, const gltexture& txt, int unit);
void set_gluniform_texture(glprogram& prog, const char* var, const char* var_on,
    const gltexture& txt, int unit);

int get_glvertexattrib_location(const glprogram& prog, const char* name);

void set_glvertexattrib(int loc, const glarraybuffer& buf, float val);
void set_glvertexattrib(int loc, const glarraybuffer& buf, const vec2f& val);
void set_glvertexattrib(int loc, const glarraybuffer& buf, const vec3f& val);
void set_glvertexattrib(int loc, const glarraybuffer& buf, const vec4f& val);

template <typename T>
inline void set_glvertexattrib(const glprogram& prog, const char* var,
    const glarraybuffer& buf, const T& val) {
    set_glvertexattrib(get_glvertexattrib_location(prog, var), buf, val);
}

void draw_glpoints(const glelementbuffer& buf, int num);
void draw_gllines(const glelementbuffer& buf, int num);
void draw_gltriangles(const glelementbuffer& buf, int num);

void draw_glimage(const gltexture& txt, vec2i imsize, vec2i winsize,
    vec2f imcenter, float imscale);

struct glwindow;

glwindow* make_glwindow(const vec2i& size, const char* title,
    void* user_pointer, function<void(glwindow*)> refresh_cb);
void      delete_glwindow(glwindow* win);

void set_drop_callback(glwindow*                               win,
    function<void(glwindow* win, const vector<string>& paths)> drop_cb);

void* get_user_pointer(glwindow* win);

vec2i get_glframebuffer_size(glwindow* win);
vec2i get_glwindow_size(glwindow* win);

bool should_glwindow_close(glwindow* win);

vec2f get_glmouse_pos(glwindow* win);
bool  get_glmouse_left(glwindow* win);
bool  get_glmouse_right(glwindow* win);
bool  get_glalt_key(glwindow* win);
bool  get_glshift_key(glwindow* win);

void process_glevents(glwindow* win, bool wait = false);
void swap_glbuffers(glwindow* win);

void init_glwidgets(glwindow* win);
bool get_glwidgets_active(glwindow* win);

void begin_glwidgets_frame(glwindow* win);
void end_glwidgets_frame(glwindow* win);

bool begin_glwidgets_window(glwindow* win, const char* title);

bool begin_header_glwidget(glwindow* win, const char* title);
void end_header_glwidget(glwindow* win);

void draw_label_glwidgets(glwindow* win, const char* lbl, const string& txt);
void draw_label_glwidgets(glwindow* win, const char* lbl, const char* fmt, ...);

bool begin_header_widget(glwindow* win, const char* label);
void end_header_widget(glwindow* win);

void draw_separator_glwidget(glwindow* win);
void continue_glwidgets_line(glwindow* win);

bool draw_button_glwidget(glwindow* win, const char* lbl);

bool draw_textinput_glwidget(glwindow* win, const char* lbl, string& val);
bool draw_slider_glwidget(
    glwindow* win, const char* lbl, float& val, float min, float max);
bool draw_slider_glwidget(
    glwindow* win, const char* lbl, vec2f& val, float min, float max);
bool draw_slider_glwidget(
    glwindow* win, const char* lbl, vec3f& val, float min, float max);
bool draw_slider_glwidget(
    glwindow* win, const char* lbl, vec4f& val, float min, float max);

bool draw_slider_glwidget(
    glwindow* win, const char* lbl, int& val, int min, int max);
bool draw_slider_glwidget(
    glwindow* win, const char* lbl, vec2i& val, int min, int max);
bool draw_slider_glwidget(
    glwindow* win, const char* lbl, vec3i& val, int min, int max);
bool draw_slider_glwidget(
    glwindow* win, const char* lbl, vec4i& val, int min, int max);

bool draw_dragger_glwidget(glwindow* win, const char* lbl, float& val,
    float speed = 1.0f, float min = 0.0f, float max = 0.0f);
bool draw_dragger_glwidget(glwindow* win, const char* lbl, vec2f& val,
    float speed = 1.0f, float min = 0.0f, float max = 0.0f);
bool draw_dragger_glwidget(glwindow* win, const char* lbl, vec3f& val,
    float speed = 1.0f, float min = 0.0f, float max = 0.0f);
bool draw_dragger_glwidget(glwindow* win, const char* lbl, vec4f& val,
    float speed = 1.0f, float min = 0.0f, float max = 0.0f);

bool draw_dragger_glwidget(glwindow* win, const char* lbl, int& val,
    float speed = 1, int min = 0, int max = 0);
bool draw_dragger_glwidget(glwindow* win, const char* lbl, vec2i& val,
    float speed = 1, int min = 0, int max = 0);
bool draw_dragger_glwidget(glwindow* win, const char* lbl, vec3i& val,
    float speed = 1, int min = 0, int max = 0);
bool draw_dragger_glwidget(glwindow* win, const char* lbl, vec4i& val,
    float speed = 1, int min = 0, int max = 0);

bool draw_checkbox_glwidget(glwindow* win, const char* lbl, bool& val);

bool draw_coloredit_glwidget(glwindow* win, const char* lbl, vec3f& val);
bool draw_coloredit_glwidget(glwindow* win, const char* lbl, vec4f& val);

bool begin_treenode_glwidget(glwindow* win, const char* lbl);
void end_treenode_glwidget(glwindow* win);

bool begin_selectabletreenode_glwidget(
    glwindow* win, const char* lbl, void*& selection, void* content);
void begin_selectabletreeleaf_glwidget(
    glwindow* win, const char* lbl, void*& selection, void* content);

bool draw_combobox_glwidget(
    glwindow* win, const char* lbl, int& idx, const vector<string>& labels);
bool draw_combobox_glwidget(
    glwindow* win, const char* lbl, string& val, const vector<string>& labels);
bool draw_combobox_glwidget(glwindow* win, const char* lbl, int& idx,
    const vector<void*>& vals, const char* (*label)(void*));
bool draw_combobox_glwidget(glwindow* win, const char* lbl, void*& val,
    const vector<void*>& vals, const char* (*label)(void*), bool include_null);

template <typename T>
bool draw_combobox_glwidget(
    glwindow* win, const char* lbl, int& idx, const vector<T*>& vals) {
    return draw_combobox_glwidget(win, lbl, idx, (const vector<void*>&)vals,
        [](void* val) { return ((T*)val)->name.c_str(); });
}

template <typename T>
bool draw_combobox_glwidget(glwindow* win, const char* lbl, T*& val,
    const vector<T*>& vals, bool include_null) {
    return draw_combobox_glwidget(win, lbl, (void*&)val,
        (const vector<void*>&)vals,
        [](void* val) { return ((T*)val)->name.c_str(); }, include_null);
}

void begin_child_glwidget(glwindow* win, const char* lbl, const vec2i& size);
void end_child_glwidget(glwindow* win);

}  // namespace ygl

#endif
