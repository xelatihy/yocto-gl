//
// Utilities to use OpenGL 3, GLFW and ImGui.
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

#include "yocto_window.h"

#include <yocto/yocto_commonio.h>

#include <array>
#include <chrono>

#include "ext/glad/glad.h"

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif
#include <GLFW/glfw3.h>

#ifdef _WIN32
#undef near
#undef far
#endif

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// UI WINDOW
// -----------------------------------------------------------------------------
namespace yocto {

inline void update_button_from_input(gui_button& button, bool pressing) {
  if (pressing) {
    // TODO(giacomo): solve this assert
    assert(button.state != gui_button::state::down);
    button.state = gui_button::state::pressing;
  } else {
    button.state = gui_button::state::releasing;
  }
}

inline void update_button_for_next_frame(gui_button& button) {
  if (button.state == gui_button::state::pressing) {
    button.state = gui_button::state::down;
  } else if (button.state == gui_button::state::releasing) {
    button.state = gui_button::state::up;
  }
}

void init_window(gui_window* win, const vec2i& size, const string& title,
    bool widgets, int widgets_width, bool widgets_left) {
  // set error callback
  glfwSetErrorCallback([](int code, const char* description) {
    printf("glfw error code: %d\n", code);
    printf("glfw error description: %s\n", description);
  });

  // init glfw
  if (!glfwInit())
    throw std::runtime_error("cannot initialize windowing system");

  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#if __APPLE__
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

  // create window
  win->title = title;
  win->win = glfwCreateWindow(size.x, size.y, title.c_str(), nullptr, nullptr);
  if (win->win == nullptr)
    throw std::runtime_error{"cannot initialize windowing system"};
  glfwMakeContextCurrent(win->win);
  glfwSwapInterval(1);  // Enable vsync

  // set user data
  glfwSetWindowUserPointer(win->win, win);

  glfwSetDropCallback(
      win->win, [](GLFWwindow* glfw, int num, const char** paths) {
        auto win = (gui_window*)glfwGetWindowUserPointer(glfw);
        win->input.dropped.resize(num);
        for (auto i = 0; i < num; i++) win->input.dropped[i] = string(paths[i]);
      });
  glfwSetKeyCallback(win->win,
      [](GLFWwindow* glfw, int key, int scancode, int action, int mods) {
        auto win   = (gui_window*)glfwGetWindowUserPointer(glfw);
        auto press = (action == GLFW_PRESS);
        update_button_from_input(win->input.key_buttons[key], press);
      });
  glfwSetCharCallback(win->win, [](GLFWwindow* glfw, unsigned int key) {
    auto win = (gui_window*)glfwGetWindowUserPointer(glfw);
    update_button_from_input(win->input.key_buttons[key], true);
  });
  glfwSetMouseButtonCallback(
      win->win, [](GLFWwindow* glfw, int button, int action, int mods) {
        auto win   = (gui_window*)glfwGetWindowUserPointer(glfw);
        auto press = (action == GLFW_PRESS);
        if (button == GLFW_MOUSE_BUTTON_LEFT) {
          update_button_from_input(win->input.mouse_left, press);
        } else if (button == GLFW_MOUSE_BUTTON_RIGHT) {
          update_button_from_input(win->input.mouse_right, press);
        } else if (button == GLFW_MOUSE_BUTTON_MIDDLE) {
          update_button_from_input(win->input.mouse_middle, press);
        }
      });
  glfwSetScrollCallback(
      win->win, [](GLFWwindow* glfw, double xoffset, double yoffset) {
        auto win          = (gui_window*)glfwGetWindowUserPointer(glfw);
        win->input.scroll = {(float)xoffset, (float)yoffset};
      });
  glfwSetWindowSizeCallback(
      win->win, [](GLFWwindow* glfw, int width, int height) {
        auto win = (gui_window*)glfwGetWindowUserPointer(glfw);
        glfwGetWindowSize(
            win->win, &win->input.window_size.x, &win->input.window_size.y);
        if (win->widgets_width) win->input.window_size.x -= win->widgets_width;
        glfwGetFramebufferSize(win->win, &win->input.framebuffer_viewport.z,
            &win->input.framebuffer_viewport.w);
        win->input.framebuffer_viewport.x = 0;
        win->input.framebuffer_viewport.y = 0;
        if (win->widgets_width) {
          auto win_size = zero2i;
          glfwGetWindowSize(win->win, &win_size.x, &win_size.y);
          auto offset = (int)(win->widgets_width *
                              (float)win->input.framebuffer_viewport.z /
                              win_size.x);
          win->input.framebuffer_viewport.z -= offset;
          if (win->widgets_left) win->input.framebuffer_viewport.x += offset;
        }
      });

  // init gl extensions
  if (!gladLoadGL())
    throw std::runtime_error{"cannot initialize OpenGL extensions"};

  // TODO(giacomo): eventually remove this stuff
  win->widgets       = true;
  win->widgets_width = widgets_width;
  win->widgets_left  = widgets_left;
}

void clear_window(gui_window* win) {
  glfwDestroyWindow(win->win);
  glfwTerminate();
  win->win = nullptr;
}

static void poll_input(gui_input& input, const gui_window* win) {
  // Clear input for next frame.
  update_button_for_next_frame(input.mouse_left);
  update_button_for_next_frame(input.mouse_right);
  for (auto& key : input.key_buttons) {
    update_button_for_next_frame(key);
  }
  input.scroll     = zero2f;
  input.mouse_last = input.mouse_pos;
  input.dropped.clear();

  // Poll new inputs
  glfwPollEvents();

  input.mouse_last = input.mouse_pos;
  auto mouse_posx = 0.0, mouse_posy = 0.0;
  glfwGetCursorPos(win->win, &mouse_posx, &mouse_posy);
  input.mouse_pos = vec2f{(float)mouse_posx, (float)mouse_posy};
  if (win->widgets_width && win->widgets_left) {
    input.mouse_pos.x -= win->widgets_width;
  }
  input.modifier_alt = glfwGetKey(win->win, GLFW_KEY_LEFT_ALT) == GLFW_PRESS ||
                       glfwGetKey(win->win, GLFW_KEY_RIGHT_ALT) == GLFW_PRESS;
  input.modifier_shift =
      glfwGetKey(win->win, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
      glfwGetKey(win->win, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS;
  input.modifier_ctrl =
      glfwGetKey(win->win, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS ||
      glfwGetKey(win->win, GLFW_KEY_RIGHT_CONTROL) == GLFW_PRESS;
  glfwGetWindowSize(win->win, &input.window_size.x, &input.window_size.y);
  if (win->widgets_width) input.window_size.x -= win->widgets_width;
  glfwGetFramebufferSize(
      win->win, &input.framebuffer_viewport.z, &input.framebuffer_viewport.w);
  input.framebuffer_viewport.x = 0;
  input.framebuffer_viewport.y = 0;
  if (win->widgets_width) {
    auto win_size = zero2i;
    glfwGetWindowSize(win->win, &win_size.x, &win_size.y);
    auto offset = (int)(win->widgets_width *
                        (float)input.framebuffer_viewport.z / win_size.x);
    input.framebuffer_viewport.z -= offset;
    if (win->widgets_left) input.framebuffer_viewport.x += offset;
  }

  // time
  input.clock_last = input.clock_now;
  input.clock_now =
      std::chrono::high_resolution_clock::now().time_since_epoch().count();
  input.time_now   = (double)input.clock_now / 1000000000.0;
  input.time_delta = (double)(input.clock_now - input.clock_last) /
                     1000000000.0;
}

void run_ui(gui_window* win, update_callback update) {
  // loop
  while (!glfwWindowShouldClose(win->win)) {
    // poll input
    poll_input(win->input, win);

    // clear framebuffer
    glClearColor(win->background.x, win->background.y, win->background.z,
        win->background.w);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // call user update function
    update(win->input, win->user_data);

    glfwSwapBuffers(win->win);
  }

  // clear
  clear_window(win);
}

void set_close(gui_window* win, bool close) {
  glfwSetWindowShouldClose(win->win, close ? GLFW_TRUE : GLFW_FALSE);
}

}  // namespace yocto
