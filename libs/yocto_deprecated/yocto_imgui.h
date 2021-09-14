//
// Yocto/ImGui: Utilities for writing native GUIs with Dear ImGui and GLFW.
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

#ifndef _YOCTO_IMGUI_
#define _YOCTO_IMGUI_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <yocto/yocto_math.h>

#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

// forward declaration
struct GLFWwindow;

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::function;
using std::string;
using std::unique_ptr;
using std::vector;

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

enum struct gui_param_type {
  // clang-format off
  value1f, value2f, value3f, value4f, 
  value1i, value2i, value3i, value4i, 
  value1s, value1b
  // clang-format on
};

struct gui_param {
  // constructors
  gui_param()
      : type{gui_param_type::value1f}
      , valuef{0, 0, 0, 0}
      , minmaxf{0, 0}
      , readonly{true} {}
  gui_param(float value, const vec2f& minmax = {0, 0}, bool readonly = false)
      : type{gui_param_type::value1f}
      , valuef{value, 0, 0, 0}
      , minmaxf{minmax}
      , readonly{readonly} {}
  gui_param(vec2f value, const vec2f& minmax = {0, 0}, bool readonly = false)
      : type{gui_param_type::value2f}
      , valuef{value.x, value.y, 0, 0}
      , minmaxf{minmax}
      , readonly{readonly} {}
  gui_param(vec3f value, const vec2f& minmax = {0, 0}, bool readonly = false)
      : type{gui_param_type::value3f}
      , valuef{value.x, value.y, value.z}
      , minmaxf{minmax}
      , readonly{readonly} {}
  gui_param(vec4f value, const vec2f& minmax = {0, 0}, bool readonly = false)
      : type{gui_param_type::value4f}
      , valuef{value.x, value.y, value.z, value.w}
      , minmaxf{minmax}
      , readonly{readonly} {}
  gui_param(vec3f value, bool color, bool readonly = false)
      : type{gui_param_type::value3f}
      , valuef{value.x, value.y, value.z, 1}
      , color{color}
      , readonly{readonly} {}
  gui_param(vec4f value, bool color, bool readonly = false)
      : type{gui_param_type::value4f}
      , valuef{value.x, value.y, value.z, value.w}
      , color{color}
      , readonly{readonly} {}
  gui_param(int value, const vec2i& minmax = {0, 0}, bool readonly = false)
      : type{gui_param_type::value1i}
      , valuei{value, 0, 0, 0}
      , minmaxi{minmax}
      , readonly{readonly} {}
  gui_param(vec2i value, const vec2i& minmax = {0, 0}, bool readonly = false)
      : type{gui_param_type::value2i}
      , valuei{value.x, value.y, 0, 0}
      , minmaxi{minmax}
      , readonly{readonly} {}
  gui_param(vec3i value, const vec2i& minmax = {0, 0}, bool readonly = false)
      : type{gui_param_type::value3i}
      , valuei{value.x, value.y, value.z, 0}
      , minmaxi{minmax}
      , readonly{readonly} {}
  gui_param(vec4i value, const vec2i& minmax = {0, 0}, bool readonly = false)
      : type{gui_param_type::value4i}
      , valuei{value.x, value.y, value.z, value.w}
      , minmaxi{minmax}
      , readonly{readonly} {}
  gui_param(bool value, bool readonly = false)
      : type{gui_param_type::value1b}, valueb{value}, readonly{readonly} {}
  gui_param(const string& value, bool readonly = false)
      : type{gui_param_type::value1s}, values{value}, readonly{readonly} {}
  gui_param(
      const string& value, const vector<string>& labels, bool readonly = false)
      : type{gui_param_type::value1s}
      , values{value}
      , labels{labels}
      , readonly{readonly} {}
  gui_param(int value, const vector<string>& labels, bool readonly = false)
      : type{gui_param_type::value1i}
      , valuei{value, 0, 0, 0}
      , labels{labels}
      , readonly{readonly} {}
  template <typename T, typename = std::enable_if_t<std::is_enum_v<T>>>
  gui_param(T value, const vector<string>& labels, bool readonly = false)
      : type{gui_param_type::value1i}
      , valuei{(int)value, 0, 0, 0}
      , labels{labels}
      , readonly{readonly} {}

  // conversions
  operator float() const {
    check_type(gui_param_type::value1f);
    return valuef.x;
  }
  operator vec2f() const {
    check_type(gui_param_type::value2f);
    return {valuef.x, valuef.y};
  }
  operator vec3f() const {
    check_type(gui_param_type::value3f);
    return {valuef.x, valuef.y, valuef.z};
  }
  operator vec4f() const {
    check_type(gui_param_type::value4f);
    return {valuef.x, valuef.y, valuef.z, valuef.w};
  }
  operator int() const {
    check_type(gui_param_type::value1i);
    return valuei.x;
  }
  operator vec2i() const {
    check_type(gui_param_type::value2i);
    return {valuei.x, valuei.y};
  }
  operator vec3i() const {
    check_type(gui_param_type::value3i);
    return {valuei.x, valuei.y, valuei.z};
  }
  operator vec4i() const {
    check_type(gui_param_type::value4i);
    return {valuei.x, valuei.y, valuei.z, valuei.w};
  }
  operator bool() const {
    check_type(gui_param_type::value1b);
    return valueb;
  }
  operator string() const {
    check_type(gui_param_type::value1s);
    return values;
  }
  template <typename T, typename = std::enable_if_t<std::is_enum_v<T>>>
  operator T() const {
    check_type(gui_param_type::value1i);
    return (T)valuei.x;
  }

  // type checking
  void check_type(gui_param_type type) const {
    if (type != this->type) throw std::invalid_argument{"bad gui type"};
  }

  // value
  gui_param_type type   = gui_param_type::value1f;
  vec4f          valuef = {0, 0, 0, 0};
  vec4i          valuei = {0, 0, 0, 0};
  bool           valueb = false;
  string         values = "";

  // display properties
  vec2f          minmaxf  = {0, 0};
  vec2i          minmaxi  = {0, 0};
  bool           color    = false;
  vector<string> labels   = {};
  bool           readonly = false;
};

struct gui_params {
  using container      = vector<pair<string, gui_param>>;
  using iterator       = container::iterator;
  using const_iterator = container::const_iterator;

  gui_params() {}

  bool   empty() const { return items.empty(); }
  size_t size() const { return items.size(); }

  gui_param& operator[](const string& key) {
    auto item = find(key);
    if (item == end()) return items.emplace_back(key, gui_param{}).second;
    return item->second;
  }
  const gui_param& operator[](const string& key) const { return at(key); }

  gui_param& at(const string& key) {
    auto item = find(key);
    if (item == end()) throw std::out_of_range{"key not found " + key};
    return item->second;
  }
  const gui_param& at(const string& key) const {
    auto item = find(key);
    if (item == end()) throw std::out_of_range{"key not found " + key};
    return item->second;
  }

  iterator find(const string& key) {
    for (auto iterator = items.begin(); iterator != items.end(); ++iterator) {
      if (iterator->first == key) return iterator;
    }
    return items.end();
  }
  const_iterator find(const string& key) const {
    for (auto iterator = items.begin(); iterator != items.end(); ++iterator) {
      if (iterator->first == key) return iterator;
    }
    return items.end();
  }

  iterator       begin() { return items.begin(); }
  iterator       end() { return items.end(); }
  const_iterator begin() const { return items.begin(); }
  const_iterator end() const { return items.end(); }

 private:
  vector<pair<string, gui_param>> items;
};

// draw params
bool draw_params(gui_window* win, const string& title, gui_params& value);

bool begin_header(gui_window* win, const char* title);
void end_header(gui_window* win);

void draw_label(gui_window* win, const char* lbl, const string& text);
void draw_label(gui_window* win, const char* lbl, int value);
void draw_label(gui_window* win, const char* lbl, bool value);

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
bool draw_coloredit(gui_window* win, const char* lbl, vec4b& value);

bool draw_hdrcoloredit(gui_window* win, const char* lbl, vec3f& value);
bool draw_hdrcoloredit(gui_window* win, const char* lbl, vec4f& value);

bool draw_combobox(gui_window* win, const char* lbl, int& idx,
    const vector<string>& labels, bool include_null = false);
bool draw_combobox(gui_window* win, const char* lbl, string& value,
    const vector<string>& labels, bool include_null = false);
bool draw_combobox(gui_window* win, const char* lbl, int& idx, int num,
    const function<string(int)>& labels, bool include_null = false);

template <typename T>
inline bool draw_combobox(gui_window* win, const char* lbl, int& idx,
    const vector<T*>& vals, bool include_null = false) {
  return draw_combobox(
      win, lbl, idx, (int)vals.size(), [&](int idx) { return vals[idx]->name; },
      include_null);
}

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
    const vector<unique_ptr<T>>& vals, bool include_null = false) {
  auto idx = -1;
  for (auto pos = 0; pos < vals.size(); pos++)
    if (vals[pos].get() == value) idx = pos;
  auto edited = draw_combobox(
      win, lbl, idx, (int)vals.size(), [&](int idx) { return vals[idx]->name; },
      include_null);
  if (edited) {
    value = idx >= 0 ? vals[idx].get() : nullptr;
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

bool draw_filedialog(gui_window* win, const char* lbl, string& path, bool save,
    const string& dirname, const string& filename, const string& filter);
bool draw_filedialog_button(gui_window* win, const char* button_lbl,
    bool button_active, const char* lbl, string& path, bool save,
    const string& dirname, const string& filename, const string& filter);

void log_info(gui_window* win, const string& msg);
void log_error(gui_window* win, const string& msg);
void clear_log(gui_window* win);
void draw_log(gui_window* win);

// Convenience funcion for params
bool draw_dragger(gui_window* win, const char* lbl, array<float, 2>& value,
    float speed = 1.0f, float min = 0.0f, float max = 0.0f);
bool draw_dragger(gui_window* win, const char* lbl, array<float, 3>& value,
    float speed = 1.0f, float min = 0.0f, float max = 0.0f);
bool draw_dragger(gui_window* win, const char* lbl, array<float, 4>& value,
    float speed = 1.0f, float min = 0.0f, float max = 0.0f);

bool draw_dragger(gui_window* win, const char* lbl, array<int, 2>& value,
    float speed = 1, int min = 0, int max = 0);
bool draw_dragger(gui_window* win, const char* lbl, array<int, 3>& value,
    float speed = 1, int min = 0, int max = 0);
bool draw_dragger(gui_window* win, const char* lbl, array<int, 4>& value,
    float speed = 1, int min = 0, int max = 0);

}  // namespace yocto

#endif
