///
/// # Yocto/Gui
///
/// A set of utilities for windows= management and immediate mode widgets
/// using GLFW and ImGui. This is just a straightfoward wrapper with the
/// main benefit of integrating ImGui in GLFW without explicit application
/// support by modying GLFW callbacks. Also provides easy to use ImGui wrappers.
///
/// This library depends in yocto_math.h
/// The library depends on GLEW for OpenGL functions on Windows and Linux,
/// GLFW for window management and dear ImGui for widgets.
///
///
/// ## History
///
/// - v 0.2: added widgets explicit IDs
/// - v 0.1: added more widgets
/// - v 0.0: initial release split from yocto_glu
///
namespace ygui {}

//
// LICENSE:
//
// Copyright (c) 2016 -- 2017 Fabio Pellacini
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

#ifndef _YUI_H_
#define _YUI_H_

#include <array>
#include <string>
#include <vector>

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// INTERFACE
// -----------------------------------------------------------------------------

///
/// Window Management and Immediate Mode Widgets
///
namespace ygui {

///
/// Forward declaration of window
///
struct window;

///
/// Text callback
///
typedef void (*text_callback)(window*, unsigned int);

///
/// Mouse callback
///
typedef void (*mouse_callback)(window*, int button, bool press, int mods);

///
/// Window refresh callback
///
typedef void (*refresh_callback)(window*);

///
/// initialize glfw
///
window* init_window(int width, int height, const std::string& title,
    void* user_pointer = nullptr);

///
/// Clear glfw
///
void clear_window(window* window);

///
/// initialize glfw
///
void set_callbacks(window* win, text_callback text_cb, mouse_callback mouse_cb,
    refresh_callback refresh_cb);

///
/// Set window title
///
void set_window_title(window* win, const std::string& title);

///
/// Wait events
///
void wait_events(window* win);

///
/// Poll events
///
void poll_events(window* win);

///
/// Swap buffers
///
void swap_buffers(window* win);

///
/// Should close
///
bool should_close(window* win);

///
/// User pointer
///
void* get_user_pointer(window* window);

///
/// Mouse button
///
int get_mouse_button(window* window);

///
/// Mouse position
///
ym::vec2i get_mouse_posi(window* window);

///
/// Mouse position
///
ym::vec2f get_mouse_posf(window* window);

///
/// Window size
///
ym::vec2i get_window_size(window* window);

///
/// Framebuffer size
///
ym::vec2i get_framebuffer_size(window* window);

///
/// Read pixels
///
std::vector<ym::vec4b> get_screenshot(
    window* win, ym::vec2i& wh, bool flipy = true, bool back = false);

///
/// Init ui
///
void init_widgets(window* win);

///
/// Clear ui
///
void clear_widgets(window* win);

///
/// Begin draw widget
///
bool begin_widgets(window* win, const std::string& title);

///
/// End draw widget
///
void end_widgets(window* win);

///
/// Whether widget are active
///
bool get_widget_active(window* win);

///
/// Horizontal separator
///
void separator_widget(window* win);

///
/// Indent widget
///
void indent_begin_widgets(window* win);

///
/// Indent widget
///
void indent_end_widgets(window* win);

///
/// Label widget
///
void label_widget(window* win, const std::string& lbl, const std::string& mesg);

///
/// Label widget
///
void label_widget(
    window* win, const std::string& lbl, int val, const char* fmt = "%d");

///
/// Label widget
///
void label_widget(window* win, const std::string& lbl, ym::vec2i val,
    const char* fmt = "[ %d, %d ]");

///
/// Label widget
///
void label_widget(window* win, const std::string& lbl, ym::vec3i val,
    const char* fmt = "[ %d, %d, %d ]");

///
/// Label widget
///
void label_widget(window* win, const std::string& lbl, ym::vec4i val,
    const char* fmt = "[ %d, %d, %d, %d ]");

///
/// Label widget
///
void label_widget(
    window* win, const std::string& lbl, float val, const char* fmt = "%f");

///
/// Label widget
///
void label_widget(window* win, const std::string& lbl, ym::vec2f val,
    const char* fmt = "[ %f, %f ]");

///
/// Label widget
///
void label_widget(window* win, const std::string& lbl, ym::vec3f val,
    const char* fmt = "[ %f, %f, %f ]");

///
/// Label widget
///
void label_widget(window* win, const std::string& lbl, ym::vec4f val,
    const char* fmt = "[ %f, %f, %f ]");

///
/// Text widget
///
bool text_widget(window* win, const std::string& lbl, char* buf, int buf_size);

///
/// Slider widget
///
bool slider_widget(window* win, const std::string& lbl, int* val, int min,
    int max, int incr = 1);

///
/// Slider widget
///
bool slider_widget(window* win, const std::string& lbl, ym::vec2i* val, int min,
    int max, int incr = 1);

///
/// Slider widget
///
bool slider_widget(window* win, const std::string& lbl, ym::vec3i* val, int min,
    int max, int incr = 1);

///
/// Slider widget
///
bool slider_widget(window* win, const std::string& lbl, ym::vec4i* val, int min,
    int max, int incr = 1);

///
/// Slider widget
///
bool slider_widget(window* win, const std::string& lbl, float* val, float min,
    float max, float incr = 1.0f);

///
/// Slider widget
///
bool slider_widget(window* win, const std::string& lbl, ym::vec2f* val,
    float min, float max, float incr = 1.0f);
///
/// Slider widget
///
bool slider_widget(window* win, const std::string& lbl, ym::vec3f* val,
    float min, float max, float incr = 1.0f);
///
/// Slider widget
///
bool slider_widget(window* win, const std::string& lbl, ym::vec4f* val,
    float min, float max, float incr = 1.0f);
///
/// Bool widget
///
bool checkbox_widget(window* win, const std::string& lbl, bool* val);

///
/// Enum widget
///
bool combo_widget(window* win, const std::string& lbl, int* val,
    const std::vector<std::pair<std::string, int>>& labels);

///
/// Enum widget
///
bool combo_widget(window* win, const std::string& lbl, void** val,
    const std::vector<std::pair<std::string, void*>>& labels);

///
/// Enum widget (assume we can convert T to int)
///
template <typename T>
inline bool combo_widget(window* win, const std::string& lbl, T* val,
    const std::vector<std::pair<std::string, T>>& labels);

///
/// Enum widget
///
template <typename T>
inline bool combo_widget(window* win, const std::string& lbl, T** val,
    const std::vector<std::pair<std::string, T*>>& labels);

///
/// List widget
///
bool list_widget(window* win, const std::string& lbl, int* val,
    const std::vector<std::pair<std::string, int>>& labels);

///
/// List widget
///
bool list_widget(window* win, const std::string& lbl, void** val,
    const std::vector<std::pair<std::string, void*>>& labels);

///
/// List widget
///
template <typename T>
inline bool list_widget(window* win, const std::string& lbl, T** val,
    const std::vector<std::pair<std::string, T*>>& labels);

///
/// Button widget
///
bool button_widget(window* win, const std::string& lbl);

///
/// Collapsible header
///
bool collapsing_header_widget(window* win, const std::string& lbl);

///
/// Start tree node
///
bool tree_begin_widget(window* win, const std::string& lbl);

///
/// End tree widget
///
void tree_end_widget(window* win);

///
/// Start selectable tree node
///
bool tree_begin_widget(
    window* win, const std::string& lbl, void** selection, void* content);

///
/// End selectable tree node
///
void tree_end_widget(window* win, void* content);

///
/// Selectable tree leaf node
///
void tree_leaf_widget(
    window* win, const std::string& lbl, void** selection, void* content);

///
/// Image widget
///
void image_widget(window* win, int tid, ym::vec2i size);

///
/// Scroll region
///
void scroll_region_begin_widget(
    window* win, const std::string& lbl, int height, bool border);

///
/// Scroll region
///
void scroll_region_end_widget(window* win);

///
/// Scroll region
///
void scroll_region_here_widget(window* win);

}  // namespace ygui

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF TEMPLATED FUNCTIONS
// -----------------------------------------------------------------------------

namespace ygui {

//
// Enum widget
//
template <typename T>
inline bool combo_widget(window* win, const std::string& lbl, T* val,
    const std::vector<std::pair<std::string, T>>& labels) {
    return combo_widget(win, lbl, (int*)val,
        (const std::vector<std::pair<std::string, int>>&)labels);
}

//
// Enum widget
//
template <typename T>
inline bool combo_widget(window* win, const std::string& lbl, T** val,
    const std::vector<std::pair<std::string, T*>>& labels) {
    return combo_widget(win, lbl, (void**)val,
        (const std::vector<std::pair<std::string, void*>>&)labels);
}

//
// List widget
//
template <typename T>
inline bool list_widget(window* win, const std::string& lbl, T** val,
    const std::vector<std::pair<std::string, T*>>& labels) {
    return list_widget(win, lbl, (void**)val,
        (const std::vector<std::pair<std::string, void*>>&)labels);
}

}  // namespace ygui

#endif
