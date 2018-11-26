//
// Yocto/OpenGL: Utilities to use OpenGL 3, GLFW and ImGui.
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

#ifndef _YOCTO_OPENGL_
#define _YOCTO_OPENGL_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "../yocto/yocto_scene.h"
#include "../yocto/yocto_utils.h"

#include <functional>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

using std::function;

}

// forward declaration
struct GLFWwindow;

// -----------------------------------------------------------------------------
// OPENGL FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

void clear_opengl_lframebuffer(const vec4f& color, bool clear_depth = true);

void set_opengl_viewport(int x, int y, int w, int h);
void set_opengl_viewport(const vec2i& size);

void set_opengl_wireframe(bool enabled);
void set_opengl_blending(bool enabled);

struct opengl_program {
    uint program_id             = 0;
    uint vertex_shader_id       = 0;
    uint fragment_shader_id     = 0;
    uint vertex_array_object_id = 0;

    operator bool() const { return (bool)program_id; }
};

bool init_opengl_program(
    opengl_program& program, const char* vertex, const char* fragment);

void delete_opengl_program(opengl_program& program);

void bind_opengl_program(opengl_program& program);
void unbind_opengl_program();

struct opengl_texture {
    uint texture_id = 0;
    int  width      = 0;
    int  height     = 0;

    operator bool() const { return (bool)texture_id; }
};

bool init_opengl_texture(opengl_texture& texture, int width, int height,
    bool as_float, bool as_srgb, bool linear, bool mipmap);

void update_opengl_texture(
    opengl_texture& texture, const image4f& img, bool mipmap);
void update_opengl_texture_region(opengl_texture& texture, const image4f& img,
    const image_region& region, bool mipmap);

inline bool init_opengl_texture(opengl_texture& texture, const image4f& img,
    bool as_float, bool linear, bool mipmap) {
    if (!init_opengl_texture(
            texture, img.width, img.height, as_float, false, linear, mipmap))
        return false;
    update_opengl_texture(texture, img, mipmap);
    return true;
}

bool init_opengl_texture(opengl_texture& texture, const image4b& img,
    bool as_srgb, bool linear, bool mipmap);
void update_opengl_texture(
    opengl_texture& texture, const image4b& img, bool mipmap);
void update_opengl_texture_region(opengl_texture& texture, const image4b& img,
    const image_region& region, bool mipmap);

inline bool init_opengl_texture(opengl_texture& texture, const image4b& img,
    bool as_srgb, bool linear, bool mipmap) {
    if (!init_opengl_texture(
            texture, img.width, img.height, false, as_srgb, linear, mipmap))
        return false;
    update_opengl_texture(texture, img, mipmap);
    return true;
}

void delete_opengl_texture(opengl_texture& texture);

struct opengl_array_buffer {
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

bool init_opengl_array_buffer(opengl_array_buffer& buffer,
    const vector<float>& data, bool dynamic = false);
bool init_opengl_array_buffer(opengl_array_buffer& buffer,
    const vector<vec2f>& data, bool dynamic = false);
bool init_opengl_array_buffer(opengl_array_buffer& buffer,
    const vector<vec3f>& data, bool dynamic = false);
bool init_opengl_array_buffer(opengl_array_buffer& buffer,
    const vector<vec4f>& data, bool dynamic = false);

void delete_opengl_array_buffer(opengl_array_buffer& buffer);

bool init_opengl_elementbuffer(opengl_elementbuffer& buffer,
    const vector<int>& data, bool dynamic = false);
bool init_opengl_elementbuffer(opengl_elementbuffer& buffer,
    const vector<vec2i>& data, bool dynamic = false);
bool init_opengl_elementbuffer(opengl_elementbuffer& buffer,
    const vector<vec3i>& data, bool dynamic = false);

void delete_opengl_elementbuffer(opengl_elementbuffer& buffer);

int get_opengl_uniform_location(
    const opengl_program& program, const char* name);

void set_opengl_uniform(int locatiom, int value);
void set_opengl_uniform(int locatiom, const vec2i& value);
void set_opengl_uniform(int locatiom, const vec3i& value);
void set_opengl_uniform(int locatiom, const vec4i& value);
void set_opengl_uniform(int locatiom, float value);
void set_opengl_uniform(int locatiom, const vec2f& value);
void set_opengl_uniform(int locatiom, const vec3f& value);
void set_opengl_uniform(int locatiom, const vec4f& value);
void set_opengl_uniform(int locatiom, const mat4f& value);
void set_opengl_uniform(int locatiom, const frame3f& value);

template <typename T>
inline void set_opengl_uniform(
    const opengl_program& program, const char* name, const T& value) {
    set_opengl_uniform(get_opengl_uniform_location(program, name), value);
}

void set_opengl_uniform_texture(
    int locatiom, const opengl_texture& texture, int unit);
void set_opengl_uniform_texture(opengl_program& program, const char* name,
    const opengl_texture& texture, int unit);
void set_opengl_uniform_texture(
    int locatiom, int locatiom_on, const opengl_texture& texture, int unit);
void set_opengl_uniform_texture(opengl_program& program, const char* name,
    const char* name_on, const opengl_texture& texture, int unit);

int get_opengl_vertexattrib_location(
    const opengl_program& program, const char* name);

void set_opengl_vertexattrib(
    int locatiom, const opengl_array_buffer& buffer, float value);
void set_opengl_vertexattrib(
    int locatiom, const opengl_array_buffer& buffer, const vec2f& value);
void set_opengl_vertexattrib(
    int locatiom, const opengl_array_buffer& buffer, const vec3f& value);
void set_opengl_vertexattrib(
    int locatiom, const opengl_array_buffer& buffer, const vec4f& value);

template <typename T>
inline void set_opengl_vertexattrib(const opengl_program& program,
    const char* name, const opengl_array_buffer& buffer, const T& value) {
    set_opengl_vertexattrib(
        get_opengl_vertexattrib_location(program, name), buffer, value);
}

void draw_opengl_points(const opengl_elementbuffer& buffer, int num);
void draw_opengl_lines(const opengl_elementbuffer& buffer, int num);
void draw_opengl_triangles(const opengl_elementbuffer& buffer, int num);

void draw_glimage(const opengl_texture& texture, int win_width, int win_height,
    const vec2f& image_center, float image_scale);
void draw_glimage_background(const opengl_texture& texture, int win_width,
    int win_height, const vec2f& image_center, float image_scale,
    float border_size = 2);

struct opengl_window;
using refresh_opengl_callback = function<void(const opengl_window&)>;
using drop_opengl_callback =
    function<void(const opengl_window&, const vector<string>&)>;

struct opengl_window {
    GLFWwindow*             win        = nullptr;
    void*                   user_ptr   = nullptr;
    refresh_opengl_callback refresh_cb = {};
    drop_opengl_callback    drop_cb    = {};
};

bool init_opengl_window(opengl_window& win, const vec2i& size,
    const string& title, void* user_pointer,
    refresh_opengl_callback refresh_cb);
void delete_opengl_window(opengl_window& win);

void set_drop_opengl_callback(opengl_window& win, drop_opengl_callback drop_cb);

void* get_opengl_user_pointer(const opengl_window& win);

vec2i get_opengl_framebuffer_size(const opengl_window& win);
vec2i get_opengl_window_size(const opengl_window& win);

bool should_opengl_window_close(const opengl_window& win);

vec2f get_opengl_mouse_pos(const opengl_window& win);
bool  get_opengl_mouse_left(const opengl_window& win);
bool  get_opengl_mouse_right(const opengl_window& win);
bool  get_opengl_alt_key(const opengl_window& win);
bool  get_opengl_shift_key(const opengl_window& win);

void process_opengl_events(const opengl_window& win, bool wait = false);
void swap_opengl_buffers(const opengl_window& win);

void init_opengl_widgets(const opengl_window& win);
bool get_opengl_widgets_active(const opengl_window& win);

void begin_opengl_widgets_frame(const opengl_window& win);
void end_opengl_widgets_frame(const opengl_window& win);

bool begin_opengl_widgets_window(const opengl_window& win, const char* title);

bool begin_header_opengl_widget(const opengl_window& win, const char* title);
void end_header_opengl_widget(const opengl_window& win);

void draw_label_opengl_widget(
    const opengl_window& win, const char* lbl, const string& texture);
void draw_label_opengl_widget(
    const opengl_window& win, const char* lbl, const char* fmt, ...);

bool begin_header_widget(const opengl_window& win, const char* label);
void end_header_widget(const opengl_window& win);

void draw_separator_opengl_widget(const opengl_window& win);
void continue_opengl_widget_line(const opengl_window& win);

bool draw_button_opengl_widget(const opengl_window& win, const char* lbl);

bool draw_textinput_opengl_widget(
    const opengl_window& win, const char* lbl, string& value);
bool draw_slider_opengl_widget(const opengl_window& win, const char* lbl,
    float& value, float min, float max);
bool draw_slider_opengl_widget(const opengl_window& win, const char* lbl,
    vec2f& value, float min, float max);
bool draw_slider_opengl_widget(const opengl_window& win, const char* lbl,
    vec3f& value, float min, float max);
bool draw_slider_opengl_widget(const opengl_window& win, const char* lbl,
    vec4f& value, float min, float max);

bool draw_slider_opengl_widget(
    const opengl_window& win, const char* lbl, int& value, int min, int max);
bool draw_slider_opengl_widget(
    const opengl_window& win, const char* lbl, vec2i& value, int min, int max);
bool draw_slider_opengl_widget(
    const opengl_window& win, const char* lbl, vec3i& value, int min, int max);
bool draw_slider_opengl_widget(
    const opengl_window& win, const char* lbl, vec4i& value, int min, int max);

bool draw_dragger_opengl_widget(const opengl_window& win, const char* lbl,
    float& value, float speed = 1.0f, float min = 0.0f, float max = 0.0f);
bool draw_dragger_opengl_widget(const opengl_window& win, const char* lbl,
    vec2f& value, float speed = 1.0f, float min = 0.0f, float max = 0.0f);
bool draw_dragger_opengl_widget(const opengl_window& win, const char* lbl,
    vec3f& value, float speed = 1.0f, float min = 0.0f, float max = 0.0f);
bool draw_dragger_opengl_widget(const opengl_window& win, const char* lbl,
    vec4f& value, float speed = 1.0f, float min = 0.0f, float max = 0.0f);

bool draw_dragger_opengl_widget(const opengl_window& win, const char* lbl,
    int& value, float speed = 1, int min = 0, int max = 0);
bool draw_dragger_opengl_widget(const opengl_window& win, const char* lbl,
    vec2i& value, float speed = 1, int min = 0, int max = 0);
bool draw_dragger_opengl_widget(const opengl_window& win, const char* lbl,
    vec3i& value, float speed = 1, int min = 0, int max = 0);
bool draw_dragger_opengl_widget(const opengl_window& win, const char* lbl,
    vec4i& value, float speed = 1, int min = 0, int max = 0);

bool draw_checkbox_opengl_widget(
    const opengl_window& win, const char* lbl, bool& value);

bool draw_coloredit_opengl_widget(
    const opengl_window& win, const char* lbl, vec3f& value);
bool draw_coloredit_opengl_widget(
    const opengl_window& win, const char* lbl, vec4f& value);

bool begin_treenode_opengl_widget(const opengl_window& win, const char* lbl);
void end_treenode_opengl_widget(const opengl_window& win);

bool begin_selectabletreenode_opengl_widget(
    const opengl_window& win, const char* lbl, bool& selected);
void begin_selectabletreeleaf_opengl_widget(
    const opengl_window& win, const char* lbl, bool& selected);

bool draw_combobox_opengl_widget(const opengl_window& win, const char* lbl,
    int& idx, const vector<string>& labels);
bool draw_combobox_opengl_widget(const opengl_window& win, const char* lbl,
    string& value, const vector<string>& labels);
bool draw_combobox_opengl_widget(const opengl_window& win, const char* lbl,
    int& idx, int num, const function<const char*(int)>& labels,
    bool include_null = false);

template <typename T>
inline bool draw_combobox_opengl_widget(const opengl_window& win,
    const char* lbl, int& idx, const vector<T*>& vals,
    bool include_null = false) {
    return draw_combobox_opengl_widget(win, lbl, idx, (int)vals.size(),
        [&](int idx) { return vals[idx]->name.c_str(); }, include_null);
}
template <typename T>
inline bool draw_combobox_opengl_widget(const opengl_window& win,
    const char* lbl, int& idx, const vector<T>& vals,
    bool include_null = false) {
    return draw_combobox_opengl_widget(win, lbl, idx, (int)vals.size(),
        [&](int idx) { return vals[idx].name.c_str(); }, include_null);
}
template <typename T>
inline bool draw_combobox_opengl_widget(const opengl_window& win,
    const char* lbl, int& idx, const deque<T>& vals,
    bool include_null = false) {
    return draw_combobox_opengl_widget(win, lbl, idx, (int)vals.size(),
        [&](int idx) { return vals[idx].name.c_str(); }, include_null);
}

void begin_child_opengl_widget(
    const opengl_window& win, const char* lbl, const vec2i& size);
void end_child_opengl_widget(const opengl_window& win);

}  // namespace yocto

#endif
