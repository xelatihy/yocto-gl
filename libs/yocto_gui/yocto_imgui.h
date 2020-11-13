//
// Yocto/ImGui: Utilities for writing native GUIs with Dear ImGui and GLFW.
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

#ifndef _YOCTO_IMGUI_
#define _YOCTO_IMGUI_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------
#include <yocto/yocto_math.h>

#include <functional>
#include <string>
#include <vector>

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
// OPENGL WIDGETS
// -----------------------------------------------------------------------------
namespace yocto {

struct gui_window;

struct gui_widgets {
  gui_window* window = nullptr;
};

gui_widgets create_imgui(gui_window* win);

bool begin_imgui(gui_widgets* widgets, const string& name,
    const vec2i& position = {0, 0}, const vec2i& size = {320, 720});
void end_imgui(gui_widgets* widgets);

bool is_active(gui_widgets* widgets);

bool begin_header(gui_widgets* widgets, const char* title);
void end_header(gui_widgets* widgets);

void draw_label(gui_widgets* widgets, const char* lbl, const string& text);

void draw_separator(gui_widgets* widgets);
void continue_line(gui_widgets* widgets);

bool draw_button(gui_widgets* widgets, const char* lbl, bool enabled = true);

bool draw_textinput(gui_widgets* widgets, const char* lbl, string& value);

bool draw_slider(
    gui_widgets* widgets, const char* lbl, float& value, float min, float max);
bool draw_slider(
    gui_widgets* widgets, const char* lbl, vec2f& value, float min, float max);
bool draw_slider(
    gui_widgets* widgets, const char* lbl, vec3f& value, float min, float max);
bool draw_slider(
    gui_widgets* widgets, const char* lbl, vec4f& value, float min, float max);

bool draw_slider(
    gui_widgets* widgets, const char* lbl, int& value, int min, int max);
bool draw_slider(
    gui_widgets* widgets, const char* lbl, vec2i& value, int min, int max);
bool draw_slider(
    gui_widgets* widgets, const char* lbl, vec3i& value, int min, int max);
bool draw_slider(
    gui_widgets* widgets, const char* lbl, vec4i& value, int min, int max);

bool draw_dragger(gui_widgets* widgets, const char* lbl, float& value,
    float speed = 1.0f, float min = 0.0f, float max = 0.0f);
bool draw_dragger(gui_widgets* widgets, const char* lbl, vec2f& value,
    float speed = 1.0f, float min = 0.0f, float max = 0.0f);
bool draw_dragger(gui_widgets* widgets, const char* lbl, vec3f& value,
    float speed = 1.0f, float min = 0.0f, float max = 0.0f);
bool draw_dragger(gui_widgets* widgets, const char* lbl, vec4f& value,
    float speed = 1.0f, float min = 0.0f, float max = 0.0f);

bool draw_dragger(gui_widgets* widgets, const char* lbl, int& value,
    float speed = 1, int min = 0, int max = 0);
bool draw_dragger(gui_widgets* widgets, const char* lbl, vec2i& value,
    float speed = 1, int min = 0, int max = 0);
bool draw_dragger(gui_widgets* widgets, const char* lbl, vec3i& value,
    float speed = 1, int min = 0, int max = 0);
bool draw_dragger(gui_widgets* widgets, const char* lbl, vec4i& value,
    float speed = 1, int min = 0, int max = 0);

bool draw_checkbox(gui_widgets* widgets, const char* lbl, bool& value);
bool draw_checkbox(
    gui_widgets* widgets, const char* lbl, bool& value, bool invert);

bool draw_coloredit(gui_widgets* widgets, const char* lbl, vec3f& value);
bool draw_coloredit(gui_widgets* widgets, const char* lbl, vec4f& value);

bool draw_hdrcoloredit(gui_widgets* widgets, const char* lbl, vec3f& value);
bool draw_hdrcoloredit(gui_widgets* widgets, const char* lbl, vec4f& value);

bool draw_combobox(gui_widgets* widgets, const char* lbl, int& idx,
    const vector<string>& labels);
bool draw_combobox(gui_widgets* widgets, const char* lbl, string& value,
    const vector<string>& labels);
bool draw_combobox(gui_widgets* widgets, const char* lbl, int& idx, int num,
    const function<string(int)>& labels, bool include_null = false);

template <typename T>
inline bool draw_combobox(gui_widgets* widgets, const char* lbl, T*& value,
    const vector<T*>& vals, bool include_null = false) {
  auto idx = -1;
  for (auto pos = 0; pos < vals.size(); pos++)
    if (vals[pos] == value) idx = pos;
  auto edited = draw_combobox(
      widgets, lbl, idx, (int)vals.size(),
      [&](int idx) { return vals[idx]->name; }, include_null);
  if (edited) {
    value = idx >= 0 ? vals[idx] : nullptr;
  }
  return edited;
}

template <typename T>
inline bool draw_combobox(gui_widgets* widgets, const char* lbl, T*& value,
    const vector<T*>& vals, const vector<string>& labels,
    bool include_null = false) {
  auto idx = -1;
  for (auto pos = 0; pos < vals.size(); pos++)
    if (vals[pos] == value) idx = pos;
  auto edited = draw_combobox(
      widgets, lbl, idx, (int)vals.size(), [&](int idx) { return labels[idx]; },
      include_null);
  if (edited) {
    value = idx >= 0 ? vals[idx] : nullptr;
  }
  return edited;
}

void draw_progressbar(gui_widgets* widgets, const char* lbl, float fraction);
void draw_progressbar(
    gui_widgets* widgets, const char* lbl, int current, int total);

void draw_histogram(
    gui_widgets* widgets, const char* lbl, const vector<float>& values);
void draw_histogram(
    gui_widgets* widgets, const char* lbl, const vector<vec2f>& values);
void draw_histogram(
    gui_widgets* widgets, const char* lbl, const vector<vec3f>& values);
void draw_histogram(
    gui_widgets* widgets, const char* lbl, const vector<vec4f>& values);

bool draw_filedialog(gui_widgets* widgets, const char* lbl, string& path,
    bool save, const string& dirname, const string& filename,
    const string& filter);
bool draw_filedialog_button(gui_widgets* widgets, const char* button_lbl,
    bool button_active, const char* lbl, string& path, bool save,
    const string& dirname, const string& filename, const string& filter);

void log_info(gui_widgets* widgets, const string& msg);
void log_error(gui_widgets* widgets, const string& msg);
void clear_log(gui_widgets* widgets);
void draw_log(gui_widgets* widgets);

}  // namespace yocto

#endif
