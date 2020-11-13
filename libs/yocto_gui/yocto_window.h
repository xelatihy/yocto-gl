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

#ifndef _YOCTO_WINDOW_
#define _YOCTO_WINDOW_

#include <yocto/yocto_math.h>

#include <array>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------

// forward declaration
struct GLFWwindow;

// -----------------------------------------------------------------------------
// UI APPLICATION
// -----------------------------------------------------------------------------
namespace yocto {

// Forward declaration of OpenGL window
struct gui_window;

enum struct gui_key : int {
  // For printable keys, just use the constructor, like gui_key('*').
  // For letters, always use upper case, like gui_key('C').

  escape        = 256,
  enter         = 257,
  tab           = 258,
  backspace     = 259,
  insert        = 260,
  _delete       = 261,
  right         = 262,
  left          = 263,
  down          = 264,
  up            = 265,
  page_up       = 266,
  page_down     = 267,
  home          = 268,
  end           = 269,
  caps_lock     = 280,
  scroll_lock   = 281,
  num_lock      = 282,
  print_screen  = 283,
  pause         = 284,
  f1            = 290,
  f2            = 291,
  f3            = 292,
  f4            = 293,
  f5            = 294,
  f6            = 295,
  f7            = 296,
  f8            = 297,
  f9            = 298,
  f10           = 299,
  f11           = 300,
  f12           = 301,
  f13           = 302,
  f14           = 303,
  f15           = 304,
  f16           = 305,
  f17           = 306,
  f18           = 307,
  f19           = 308,
  f20           = 309,
  f21           = 310,
  f22           = 311,
  f23           = 312,
  f24           = 313,
  f25           = 314,
  kp_0          = 320,
  kp_1          = 321,
  kp_2          = 322,
  kp_3          = 323,
  kp_4          = 324,
  kp_5          = 325,
  kp_6          = 326,
  kp_7          = 327,
  kp_8          = 328,
  kp_9          = 329,
  kp_decimal    = 330,
  kp_divide     = 331,
  kp_multiply   = 332,
  kp_subtract   = 333,
  kp_add        = 334,
  kp_enter      = 335,
  kp_equal      = 336,
  left_shift    = 340,
  left_control  = 341,
  left_alt      = 342,
  left_super    = 343,
  right_shift   = 344,
  right_control = 345,
  right_alt     = 346,
  right_super   = 347,
  menu          = 348,
  world_1       = 161,  //  non-us #1
  world_2       = 162   //  non-us #2

};

struct gui_button {
  enum struct state : unsigned char {
    up = 0,
    down,
    releasing,
    pressing,
  };
  state state = state::up;

  operator bool() const { return state == state::down; }
};

// Input state
struct gui_input {
  // Ever-changing data
  vec2f    mouse_pos  = {0, 0};  // position excluding gui
  vec2f    mouse_last = {0, 0};  // last mouse position excluding gui
  uint64_t clock_now  = 0;       // clock now
  uint64_t clock_last = 0;       // clock last
  double   time_now   = 0;       // time now
  double   time_delta = 0;       // time delta
  int      frame      = 0;

  gui_button mouse_left   = {};
  gui_button mouse_middle = {};
  gui_button mouse_right  = {};
  vec2f      scroll       = {0, 0};  // scroll input

  bool modifier_alt   = false;  // alt modifier
  bool modifier_ctrl  = false;  // ctrl modifier
  bool modifier_shift = false;  // shift modifier

  vec2i window_size          = {0, 0};
  vec4i framebuffer_viewport = {0, 0, 0, 0};
  vec2i framebuffer_size     = {0, 0};

  std::vector<std::string>    dropped     = {};
  std::array<gui_button, 512> key_buttons = {};
};

}  // namespace yocto

// -----------------------------------------------------------------------------
// UI WINDOW
// -----------------------------------------------------------------------------
namespace yocto {

// OpenGL window wrapper
struct gui_window {
  GLFWwindow* win           = nullptr;
  std::string title         = "";
  bool        widgets       = false;
  int         widgets_width = 0;
  bool        widgets_left  = true;
  gui_input   input         = {};
  vec4f       background    = {0.15f, 0.15f, 0.15f, 1.0f};
  void*       user_data     = nullptr;
};

// Windows initialization
void init_window(gui_window* win, const vec2i& size, const std::string& title,
    bool widgets, int widgets_width = 320, bool widgets_left = true);

// Window cleanup
void clear_window(gui_window* win);

// Callback
typedef void (*update_callback)(const gui_input&, void*);

// Run loop
void run_ui(gui_window* win, update_callback update);

void set_close(gui_window* win, bool close);

}  // namespace yocto

#endif
