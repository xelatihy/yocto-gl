//
// Utilities to use OpenGL 3, GLFW and ImGui.
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

#include "yocto_imgui.h"

#include <yocto/yocto_color.h>

#include <algorithm>
#include <array>
#include <cstdarg>
#include <filesystem>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <unordered_map>
#include <utility>

#include "ext/glad/glad.h"

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif
#include <GLFW/glfw3.h>

#include "ext/imgui/imgui.h"
#include "ext/imgui/imgui_impl_glfw.h"
#include "ext/imgui/imgui_impl_opengl3.h"
#include "ext/imgui/imgui_internal.h"

#ifdef _WIN32
#undef near
#undef far
#endif

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::mutex;
using std::pair;
using std::unordered_map;
using namespace std::string_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// PATH UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Make a path from a utf8 string
static std::filesystem::path make_path(const string& filename) {
  return std::filesystem::u8path(filename);
}

// Normalize path
static string normalize_path(const string& filename) {
  return make_path(filename).generic_u8string();
}

// Get extension (including .)
static string path_extension(const string& filename) {
  return make_path(filename).extension().u8string();
}

// Get filename without directory.
static string path_filename(const string& filename) {
  return make_path(filename).filename().u8string();
}

// Get filename without directory and extension.
static string path_basename(const string& filename) {
  return make_path(filename).stem().u8string();
}

// Joins paths
static string path_join(const string& patha, const string& pathb) {
  return (make_path(patha) / make_path(pathb)).generic_u8string();
}

// Check if a file can be opened for reading.
static bool path_exists(const string& filename) {
  return exists(make_path(filename));
}

// Check if a file is a directory
static bool path_isdir(const string& filename) {
  return is_directory(make_path(filename));
}

// List the contents of a directory
static vector<string> list_directory(const string& filename) {
  auto entries = vector<string>{};
  for (auto entry : std::filesystem::directory_iterator(make_path(filename))) {
    entries.push_back(entry.path().generic_u8string());
  }
  return entries;
}

// Get the current directory
static string path_current() {
  return std::filesystem::current_path().u8string();
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// UI APPLICATION
// -----------------------------------------------------------------------------
namespace yocto {

// run the user interface with the give callbacks
void run_ui(const vec2i& size, const string& title,
    const gui_callbacks& callbacks, int widgets_width, bool widgets_left) {
  auto win_guard = std::make_unique<gui_window>();
  auto win       = win_guard.get();
  init_window(win, size, title, (bool)callbacks.widgets_cb, widgets_width,
      widgets_left);

  set_init_callback(win, callbacks.init_cb);
  set_clear_callback(win, callbacks.clear_cb);
  set_draw_callback(win, callbacks.draw_cb);
  set_widgets_callback(win, callbacks.widgets_cb);
  set_drop_callback(win, callbacks.drop_cb);
  set_key_callback(win, callbacks.key_cb);
  set_char_callback(win, callbacks.char_cb);
  set_click_callback(win, callbacks.click_cb);
  set_scroll_callback(win, callbacks.scroll_cb);
  set_update_callback(win, callbacks.update_cb);
  set_uiupdate_callback(win, callbacks.uiupdate_cb);

  // run ui
  run_ui(win);

  // clear
  clear_window(win);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// UI WINDOW
// -----------------------------------------------------------------------------
namespace yocto {

static void draw_window(gui_window* win) {
  glClearColor(win->background.x, win->background.y, win->background.z,
      win->background.w);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  if (win->draw_cb) win->draw_cb(win, win->input);
  if (win->widgets_cb) {
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    auto window = zero2i;
    glfwGetWindowSize(win->win, &window.x, &window.y);
    if (win->widgets_left) {
      ImGui::SetNextWindowPos({0, 0});
      ImGui::SetNextWindowSize({(float)win->widgets_width, (float)window.y});
    } else {
      ImGui::SetNextWindowPos({(float)(window.x - win->widgets_width), 0});
      ImGui::SetNextWindowSize({(float)win->widgets_width, (float)window.y});
    }
    ImGui::SetNextWindowCollapsed(false);
    ImGui::SetNextWindowBgAlpha(1);
    if (ImGui::Begin(win->title.c_str(), nullptr,
            ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
                ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse |
                ImGuiWindowFlags_NoSavedSettings)) {
      win->widgets_cb(win, win->input);
    }
    ImGui::End();
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
  }
  glfwSwapBuffers(win->win);
}

void init_window(gui_window* win, const vec2i& size, const string& title,
    bool widgets, int widgets_width, bool widgets_left) {
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

  // set callbacks
  glfwSetWindowRefreshCallback(win->win, [](GLFWwindow* glfw) {
    auto win = (gui_window*)glfwGetWindowUserPointer(glfw);
    draw_window(win);
  });
  glfwSetDropCallback(
      win->win, [](GLFWwindow* glfw, int num, const char** paths) {
        auto win = (gui_window*)glfwGetWindowUserPointer(glfw);
        if (win->drop_cb) {
          auto pathv = vector<string>();
          for (auto i = 0; i < num; i++) pathv.push_back(paths[i]);
          win->drop_cb(win, pathv, win->input);
        }
      });
  glfwSetKeyCallback(win->win,
      [](GLFWwindow* glfw, int key, int scancode, int action, int mods) {
        auto win = (gui_window*)glfwGetWindowUserPointer(glfw);
        if (win->key_cb) win->key_cb(win, key, (bool)action, win->input);
      });
  glfwSetCharCallback(win->win, [](GLFWwindow* glfw, unsigned int key) {
    auto win = (gui_window*)glfwGetWindowUserPointer(glfw);
    if (win->char_cb) win->char_cb(win, key, win->input);
  });
  glfwSetMouseButtonCallback(
      win->win, [](GLFWwindow* glfw, int button, int action, int mods) {
        auto win = (gui_window*)glfwGetWindowUserPointer(glfw);
        if (win->click_cb)
          win->click_cb(
              win, button == GLFW_MOUSE_BUTTON_LEFT, (bool)action, win->input);
      });
  glfwSetScrollCallback(
      win->win, [](GLFWwindow* glfw, double xoffset, double yoffset) {
        auto win = (gui_window*)glfwGetWindowUserPointer(glfw);
        if (win->scroll_cb) win->scroll_cb(win, (float)yoffset, win->input);
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

  // widgets
  if (widgets) {
    ImGui::CreateContext();
    ImGui::GetIO().IniFilename       = nullptr;
    ImGui::GetStyle().WindowRounding = 0;
    ImGui_ImplGlfw_InitForOpenGL(win->win, true);
#ifndef __APPLE__
    ImGui_ImplOpenGL3_Init();
#else
    ImGui_ImplOpenGL3_Init("#version 330");
#endif
    ImGui::StyleColorsDark();
    win->widgets_width = widgets_width;
    win->widgets_left  = widgets_left;
  }
}

void clear_window(gui_window* win) {
  glfwDestroyWindow(win->win);
  glfwTerminate();
  win->win = nullptr;
}

// Run loop
void run_ui(gui_window* win) {
  // init
  if (win->init_cb) win->init_cb(win, win->input);

  // loop
  while (!glfwWindowShouldClose(win->win)) {
    // update input
    win->input.mouse_last = win->input.mouse_pos;
    auto mouse_posx = 0.0, mouse_posy = 0.0;
    glfwGetCursorPos(win->win, &mouse_posx, &mouse_posy);
    win->input.mouse_pos = vec2f{(float)mouse_posx, (float)mouse_posy};
    if (win->widgets_width && win->widgets_left)
      win->input.mouse_pos.x -= win->widgets_width;
    win->input.mouse_left = glfwGetMouseButton(
                                win->win, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
    win->input.mouse_right =
        glfwGetMouseButton(win->win, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS;
    win->input.modifier_alt =
        glfwGetKey(win->win, GLFW_KEY_LEFT_ALT) == GLFW_PRESS ||
        glfwGetKey(win->win, GLFW_KEY_RIGHT_ALT) == GLFW_PRESS;
    win->input.modifier_shift =
        glfwGetKey(win->win, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
        glfwGetKey(win->win, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS;
    win->input.modifier_ctrl =
        glfwGetKey(win->win, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS ||
        glfwGetKey(win->win, GLFW_KEY_RIGHT_CONTROL) == GLFW_PRESS;
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
    if (win->widgets_width) {
      auto io                   = &ImGui::GetIO();
      win->input.widgets_active = io->WantTextInput || io->WantCaptureMouse ||
                                  io->WantCaptureKeyboard;
    }

    // time
    win->input.clock_last = win->input.clock_now;
    win->input.clock_now =
        std::chrono::high_resolution_clock::now().time_since_epoch().count();
    win->input.time_now = (double)win->input.clock_now / 1000000000.0;
    win->input.time_delta =
        (double)(win->input.clock_now - win->input.clock_last) / 1000000000.0;

    // update ui
    if (win->uiupdate_cb && !win->input.widgets_active)
      win->uiupdate_cb(win, win->input);

    // update
    if (win->update_cb) win->update_cb(win, win->input);

    // draw
    draw_window(win);

    // event hadling
    glfwPollEvents();
  }

  // clear
  if (win->clear_cb) win->clear_cb(win, win->input);
}

void set_init_callback(gui_window* win, init_callback cb) { win->init_cb = cb; }
void set_clear_callback(gui_window* win, clear_callback cb) {
  win->clear_cb = cb;
}
void set_draw_callback(gui_window* win, draw_callback cb) { win->draw_cb = cb; }
void set_widgets_callback(gui_window* win, widgets_callback cb) {
  win->widgets_cb = cb;
}
void set_drop_callback(gui_window* win, drop_callback drop_cb) {
  win->drop_cb = drop_cb;
}
void set_key_callback(gui_window* win, key_callback cb) { win->key_cb = cb; }
void set_char_callback(gui_window* win, char_callback cb) { win->char_cb = cb; }
void set_click_callback(gui_window* win, click_callback cb) {
  win->click_cb = cb;
}
void set_scroll_callback(gui_window* win, scroll_callback cb) {
  win->scroll_cb = cb;
}
void set_uiupdate_callback(gui_window* win, uiupdate_callback cb) {
  win->uiupdate_cb = cb;
}
void set_update_callback(gui_window* win, update_callback cb) {
  win->update_cb = cb;
}

void set_close(gui_window* win, bool close) {
  glfwSetWindowShouldClose(win->win, close ? GLFW_TRUE : GLFW_FALSE);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// OPENGL WIDGETS
// -----------------------------------------------------------------------------
namespace yocto {

void init_glwidgets(gui_window* win, int width, bool left) {
  // init widgets
  ImGui::CreateContext();
  ImGui::GetIO().IniFilename       = nullptr;
  ImGui::GetStyle().WindowRounding = 0;
  ImGui_ImplGlfw_InitForOpenGL(win->win, true);
#ifndef __APPLE__
  ImGui_ImplOpenGL3_Init();
#else
  ImGui_ImplOpenGL3_Init("#version 330");
#endif
  ImGui::StyleColorsDark();
  win->widgets_width = width;
  win->widgets_left  = left;
}

bool begin_header(gui_window* win, const char* lbl) {
  if (!ImGui::CollapsingHeader(lbl)) return false;
  ImGui::PushID(lbl);
  return true;
}
void end_header(gui_window* win) { ImGui::PopID(); }

void open_glmodal(gui_window* win, const char* lbl) { ImGui::OpenPopup(lbl); }
void clear_glmodal(gui_window* win) { ImGui::CloseCurrentPopup(); }
bool begin_glmodal(gui_window* win, const char* lbl) {
  return ImGui::BeginPopupModal(lbl);
}
void end_glmodal(gui_window* win) { ImGui::EndPopup(); }
bool is_glmodal_open(gui_window* win, const char* lbl) {
  return ImGui::IsPopupOpen(lbl);
}

struct filedialog_state {
  string                     dirname       = "";
  string                     filename      = "";
  vector<pair<string, bool>> entries       = {};
  bool                       save          = false;
  bool                       remove_hidden = true;
  string                     filter        = "";
  vector<string>             extensions    = {};

  filedialog_state() = default;
  filedialog_state(const string& dirname, const string& filename,
      const string& filter, bool save) {
    set(dirname, filename, filter, save);
  }

  void set(const string& dirname, const string& filename, const string& filter,
      bool save) {
    this->save = save;
    _set_filter(filter);
    _set_dirname(dirname);
    _set_filename(filename);
  }

  void _set_dirname(const string& name) {
    if (path_exists(name) && path_isdir(name)) {
      dirname = name;
    } else if (path_exists(dirname) && path_isdir(dirname)) {
      // leave it like this
    } else {
      dirname = path_current();
    }
    dirname = normalize_path(dirname);
    entries.clear();
    for (auto entry : list_directory(dirname)) {
      if (remove_hidden && path_basename(entry)[0] == '.') continue;
      if (path_isdir(entry)) {
        entries.push_back({path_filename(entry) + "/", true});
      } else {
        entries.push_back({path_filename(entry), false});
      }
    }
    std::sort(entries.begin(), entries.end(), [](auto& a, auto& b) {
      if (a.second == b.second) return a.first < b.first;
      return a.second;
    });
  }

  void _set_filename(const string& name) {
    filename = name;
    if (filename.empty()) return;
    auto ext = path_extension(filename);
    if (std::find(extensions.begin(), extensions.end(), ext) ==
        extensions.end()) {
      filename = "";
      return;
    }
    if (!save && !path_exists(path_join(dirname, filename))) {
      filename = "";
      return;
    }
  }

  void _set_filter(const string& flt) {
    auto globs = vector<string>{""};
    for (auto i = 0; i < flt.size(); i++) {
      if (flt[i] == ';') {
        globs.push_back("");
      } else {
        globs.back() += flt[i];
      }
    }
    filter = "";
    extensions.clear();
    for (auto pattern : globs) {
      if (pattern == "") continue;
      auto ext = path_extension(pattern);
      if (ext != "") {
        extensions.push_back(ext);
        filter += (filter == "") ? ("*." + ext) : (";*." + ext);
      }
    }
  }

  void select(int idx) {
    if (entries[idx].second) {
      set(path_join(dirname, entries[idx].first), filename, filter, save);
    } else {
      set(dirname, entries[idx].first, filter, save);
    }
  }

  string get_path() const { return path_join(dirname, filename); }
};

bool draw_filedialog(gui_window* win, const char* lbl, string& path, bool save,
    const string& dirname, const string& filename, const string& filter) {
  static auto states = unordered_map<string, filedialog_state>{};
  ImGui::SetNextWindowSize({500, 300}, ImGuiCond_FirstUseEver);
  if (ImGui::BeginPopupModal(lbl)) {
    if (states.find(lbl) == states.end()) {
      states[lbl] = filedialog_state{dirname, filename, filter, save};
    }
    auto& state      = states.at(lbl);
    auto  dir_buffer = array<char, 1024>{};
    snprintf(dir_buffer.data(), dir_buffer.size(), "%s", state.dirname.c_str());
    if (ImGui::InputText("dir", dir_buffer.data(), sizeof(dir_buffer))) {
      state.set(dir_buffer.data(), state.filename, state.filter, save);
    }
    auto current_item = -1;
    if (ImGui::ListBox(
            "entries", &current_item,
            [](void* data, int idx, const char** out_text) -> bool {
              auto& state = *(filedialog_state*)data;
              *out_text   = state.entries[idx].first.c_str();
              return true;
            },
            &state, (int)state.entries.size())) {
      state.select(current_item);
    }
    auto file_buffer = array<char, 1024>{};
    snprintf(
        file_buffer.data(), file_buffer.size(), "%s", state.filename.c_str());
    if (ImGui::InputText("file", file_buffer.data(), file_buffer.size())) {
      state.set(state.dirname, file_buffer.data(), state.filter, save);
    }
    auto filter_buffer = array<char, 1024>{};
    snprintf(
        filter_buffer.data(), filter_buffer.size(), "%s", state.filter.c_str());
    if (ImGui::InputText(
            "filter", filter_buffer.data(), filter_buffer.size())) {
      state.set(state.dirname, state.filename, filter_buffer.data(), save);
    }
    auto ok = false, exit = false;
    if (ImGui::Button("Ok")) {
      path = state.dirname + state.filename;
      ok   = true;
      exit = true;
    }
    ImGui::SameLine();
    if (ImGui::Button("Cancel")) {
      exit = true;
    }
    if (exit) {
      ImGui::CloseCurrentPopup();
      states.erase(lbl);
    }
    ImGui::EndPopup();
    return ok;
  } else {
    return false;
  }
}

bool draw_filedialog_button(gui_window* win, const char* button_lbl,
    bool button_active, const char* lbl, string& path, bool save,
    const string& dirname, const string& filename, const string& filter) {
  if (is_glmodal_open(win, lbl)) {
    draw_button(win, button_lbl, button_active);
    return draw_filedialog(win, lbl, path, save, dirname, filename, filter);
  } else {
    if (draw_button(win, button_lbl, button_active)) {
      open_glmodal(win, lbl);
    }
    return false;
  }
}

bool draw_button(gui_window* win, const char* lbl, bool enabled) {
  if (enabled) {
    return ImGui::Button(lbl);
  } else {
    ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
    ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
    auto ok = ImGui::Button(lbl);
    ImGui::PopItemFlag();
    ImGui::PopStyleVar();
    return ok;
  }
}

void draw_label(gui_window* win, const char* lbl, const string& label) {
  ImGui::LabelText(lbl, "%s", label.c_str());
}

void draw_separator(gui_window* win) { ImGui::Separator(); }

void continue_line(gui_window* win) { ImGui::SameLine(); }

bool draw_textinput(gui_window* win, const char* lbl, string& value) {
  auto buffer = array<char, 4096>{};
  auto num    = 0;
  for (auto c : value) buffer[num++] = c;
  buffer[num] = 0;
  auto edited = ImGui::InputText(lbl, buffer.data(), buffer.size());
  if (edited) value = buffer.data();
  return edited;
}

bool draw_slider(
    gui_window* win, const char* lbl, float& value, float min, float max) {
  return ImGui::SliderFloat(lbl, &value, min, max);
}
bool draw_slider(
    gui_window* win, const char* lbl, vec2f& value, float min, float max) {
  return ImGui::SliderFloat2(lbl, &value.x, min, max);
}
bool draw_slider(
    gui_window* win, const char* lbl, vec3f& value, float min, float max) {
  return ImGui::SliderFloat3(lbl, &value.x, min, max);
}
bool draw_slider(
    gui_window* win, const char* lbl, vec4f& value, float min, float max) {
  return ImGui::SliderFloat4(lbl, &value.x, min, max);
}

bool draw_slider(
    gui_window* win, const char* lbl, int& value, int min, int max) {
  return ImGui::SliderInt(lbl, &value, min, max);
}
bool draw_slider(
    gui_window* win, const char* lbl, vec2i& value, int min, int max) {
  return ImGui::SliderInt2(lbl, &value.x, min, max);
}
bool draw_slider(
    gui_window* win, const char* lbl, vec3i& value, int min, int max) {
  return ImGui::SliderInt3(lbl, &value.x, min, max);
}
bool draw_slider(
    gui_window* win, const char* lbl, vec4i& value, int min, int max) {
  return ImGui::SliderInt4(lbl, &value.x, min, max);
}

bool draw_dragger(gui_window* win, const char* lbl, float& value, float speed,
    float min, float max) {
  return ImGui::DragFloat(lbl, &value, speed, min, max);
}
bool draw_dragger(gui_window* win, const char* lbl, vec2f& value, float speed,
    float min, float max) {
  return ImGui::DragFloat2(lbl, &value.x, speed, min, max);
}
bool draw_dragger(gui_window* win, const char* lbl, vec3f& value, float speed,
    float min, float max) {
  return ImGui::DragFloat3(lbl, &value.x, speed, min, max);
}
bool draw_dragger(gui_window* win, const char* lbl, vec4f& value, float speed,
    float min, float max) {
  return ImGui::DragFloat4(lbl, &value.x, speed, min, max);
}

bool draw_dragger(gui_window* win, const char* lbl, int& value, float speed,
    int min, int max) {
  return ImGui::DragInt(lbl, &value, speed, min, max);
}
bool draw_dragger(gui_window* win, const char* lbl, vec2i& value, float speed,
    int min, int max) {
  return ImGui::DragInt2(lbl, &value.x, speed, min, max);
}
bool draw_dragger(gui_window* win, const char* lbl, vec3i& value, float speed,
    int min, int max) {
  return ImGui::DragInt3(lbl, &value.x, speed, min, max);
}
bool draw_dragger(gui_window* win, const char* lbl, vec4i& value, float speed,
    int min, int max) {
  return ImGui::DragInt4(lbl, &value.x, speed, min, max);
}

bool draw_dragger(gui_window* win, const char* lbl, array<float, 2>& value,
    float speed, float min, float max) {
  return ImGui::DragFloat2(lbl, value.data(), speed, min, max);
}
bool draw_dragger(gui_window* win, const char* lbl, array<float, 3>& value,
    float speed, float min, float max) {
  return ImGui::DragFloat3(lbl, value.data(), speed, min, max);
}
bool draw_dragger(gui_window* win, const char* lbl, array<float, 4>& value,
    float speed, float min, float max) {
  return ImGui::DragFloat4(lbl, value.data(), speed, min, max);
}

bool draw_dragger(gui_window* win, const char* lbl, array<int, 2>& value,
    float speed, int min, int max) {
  return ImGui::DragInt2(lbl, value.data(), speed, min, max);
}
bool draw_dragger(gui_window* win, const char* lbl, array<int, 3>& value,
    float speed, int min, int max) {
  return ImGui::DragInt3(lbl, value.data(), speed, min, max);
}
bool draw_dragger(gui_window* win, const char* lbl, array<int, 4>& value,
    float speed, int min, int max) {
  return ImGui::DragInt4(lbl, value.data(), speed, min, max);
}

bool draw_checkbox(gui_window* win, const char* lbl, bool& value) {
  return ImGui::Checkbox(lbl, &value);
}
bool draw_checkbox(gui_window* win, const char* lbl, bool& value, bool invert) {
  if (!invert) {
    return draw_checkbox(win, lbl, value);
  } else {
    auto inverted = !value;
    auto edited   = ImGui::Checkbox(lbl, &inverted);
    if (edited) value = !inverted;
    return edited;
  }
}

bool draw_coloredit(gui_window* win, const char* lbl, vec3f& value) {
  auto flags = ImGuiColorEditFlags_Float;
  return ImGui::ColorEdit3(lbl, &value.x, flags);
}

bool draw_coloredit(gui_window* win, const char* lbl, vec4f& value) {
  auto flags = ImGuiColorEditFlags_Float;
  return ImGui::ColorEdit4(lbl, &value.x, flags);
}

bool draw_hdrcoloredit(gui_window* win, const char* lbl, vec3f& value) {
  auto color    = value;
  auto exposure = 0.0f;
  auto scale    = max(color);
  if (scale > 1) {
    color /= scale;
    exposure = log2(scale);
  }
  auto edit_exposure = draw_slider(
      win, (lbl + " [exp]"s).c_str(), exposure, 0, 10);
  auto edit_color = draw_coloredit(win, (lbl + " [col]"s).c_str(), color);
  if (edit_exposure || edit_color) {
    value = color * exp2(exposure);
    return true;
  } else {
    return false;
  }
}
bool draw_hdrcoloredit(gui_window* win, const char* lbl, vec4f& value) {
  auto color    = value;
  auto exposure = 0.0f;
  auto scale    = max(xyz(color));
  if (scale > 1) {
    color.x /= scale;
    color.y /= scale;
    color.z /= scale;
    exposure = log2(scale);
  }
  auto edit_exposure = draw_slider(
      win, (lbl + " [exp]"s).c_str(), exposure, 0, 10);
  auto edit_color = draw_coloredit(win, (lbl + " [col]"s).c_str(), color);
  if (edit_exposure || edit_color) {
    value.x = color.x * exp2(exposure);
    value.y = color.y * exp2(exposure);
    value.z = color.z * exp2(exposure);
    value.w = color.w;
    return true;
  } else {
    return false;
  }
}

bool draw_coloredit(gui_window* win, const char* lbl, vec4b& value) {
  auto valuef = byte_to_float(value);
  if (ImGui::ColorEdit4(lbl, &valuef.x)) {
    value = float_to_byte(valuef);
    return true;
  } else {
    return false;
  }
}

bool draw_combobox(gui_window* win, const char* lbl, int& value,
    const vector<string>& labels, bool include_null) {
  if (!ImGui::BeginCombo(lbl, value >= 0 ? labels.at(value).c_str() : "<none>"))
    return false;
  auto old_val = value;
  if (include_null) {
    ImGui::PushID(100000);
    if (ImGui::Selectable("<none>", value < 0)) value = -1;
    if (value < 0) ImGui::SetItemDefaultFocus();
    ImGui::PopID();
  }
  for (auto i = 0; i < labels.size(); i++) {
    ImGui::PushID(i);
    if (ImGui::Selectable(labels[i].c_str(), value == i)) value = i;
    if (value == i) ImGui::SetItemDefaultFocus();
    ImGui::PopID();
  }
  ImGui::EndCombo();
  return value != old_val;
}

bool draw_combobox(gui_window* win, const char* lbl, string& value,
    const vector<string>& labels, bool include_null) {
  if (!ImGui::BeginCombo(lbl, value.c_str())) return false;
  auto old_val = value;
  if (include_null) {
    ImGui::PushID(100000);
    if (ImGui::Selectable("<none>", value.empty())) value = "";
    if (value.empty()) ImGui::SetItemDefaultFocus();
    ImGui::PopID();
  }
  for (auto i = 0; i < labels.size(); i++) {
    ImGui::PushID(i);
    if (ImGui::Selectable(labels[i].c_str(), value == labels[i]))
      value = labels[i];
    if (value == labels[i]) ImGui::SetItemDefaultFocus();
    ImGui::PopID();
  }
  ImGui::EndCombo();
  return value != old_val;
}

bool draw_combobox(gui_window* win, const char* lbl, int& idx, int num,
    const function<string(int)>& labels, bool include_null) {
  if (num <= 0) idx = -1;
  if (!ImGui::BeginCombo(lbl, idx >= 0 ? labels(idx).c_str() : "<none>"))
    return false;
  auto old_idx = idx;
  if (include_null) {
    ImGui::PushID(100000);
    if (ImGui::Selectable("<none>", idx < 0)) idx = -1;
    if (idx < 0) ImGui::SetItemDefaultFocus();
    ImGui::PopID();
  }
  for (auto i = 0; i < num; i++) {
    ImGui::PushID(i);
    if (ImGui::Selectable(labels(i).c_str(), idx == i)) idx = i;
    if (idx == i) ImGui::SetItemDefaultFocus();
    ImGui::PopID();
  }
  ImGui::EndCombo();
  return idx != old_idx;
}

void draw_progressbar(gui_window* win, const char* lbl, float fraction) {
  ImGui::PushStyleColor(ImGuiCol_PlotHistogram, ImVec4(0.5, 0.5, 1, 0.25));
  ImGui::ProgressBar(fraction, ImVec2(0.0f, 0.0f));
  ImGui::SameLine(0.0f, ImGui::GetStyle().ItemInnerSpacing.x);
  ImGui::Text(lbl, ImVec2(0.0f, 0.0f));
  ImGui::PopStyleColor(1);
}

void draw_progressbar(
    gui_window* win, const char* lbl, int current, int total) {
  auto overlay = std::to_string(current) + "/" + std::to_string(total);
  ImGui::PushStyleColor(ImGuiCol_PlotHistogram, ImVec4(0.5, 0.5, 1, 0.25));
  ImGui::ProgressBar(
      (float)current / (float)total, ImVec2(0.0f, 0.0f), overlay.c_str());
  ImGui::SameLine(0.0f, ImGui::GetStyle().ItemInnerSpacing.x);
  ImGui::Text(lbl, ImVec2(0.0f, 0.0f));
  ImGui::PopStyleColor(1);
}

void draw_histogram(
    gui_window* win, const char* lbl, const float* values, int count) {
  ImGui::PlotHistogram(lbl, values, count);
}
void draw_histogram(
    gui_window* win, const char* lbl, const vector<float>& values) {
  ImGui::PlotHistogram(lbl, values.data(), (int)values.size(), 0, nullptr,
      flt_max, flt_max, {0, 0}, 4);
}
void draw_histogram(
    gui_window* win, const char* lbl, const vector<vec2f>& values) {
  ImGui::PlotHistogram((lbl + " x"s).c_str(), (const float*)values.data() + 0,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec2f));
  ImGui::PlotHistogram((lbl + " y"s).c_str(), (const float*)values.data() + 1,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec2f));
}
void draw_histogram(
    gui_window* win, const char* lbl, const vector<vec3f>& values) {
  ImGui::PlotHistogram((lbl + " x"s).c_str(), (const float*)values.data() + 0,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec3f));
  ImGui::PlotHistogram((lbl + " y"s).c_str(), (const float*)values.data() + 1,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec3f));
  ImGui::PlotHistogram((lbl + " z"s).c_str(), (const float*)values.data() + 2,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec3f));
}
void draw_histogram(
    gui_window* win, const char* lbl, const vector<vec4f>& values) {
  ImGui::PlotHistogram((lbl + " x"s).c_str(), (const float*)values.data() + 0,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec4f));
  ImGui::PlotHistogram((lbl + " y"s).c_str(), (const float*)values.data() + 1,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec4f));
  ImGui::PlotHistogram((lbl + " z"s).c_str(), (const float*)values.data() + 2,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec4f));
  ImGui::PlotHistogram((lbl + " w"s).c_str(), (const float*)values.data() + 3,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec4f));
}

// https://github.com/ocornut/imgui/issues/300
struct ImGuiAppLog {
  ImGuiTextBuffer Buf;
  ImGuiTextFilter Filter;
  ImVector<int>   LineOffsets;  // Index to lines offset
  bool            ScrollToBottom;

  void Clear() {
    Buf.clear();
    LineOffsets.clear();
  }

  void AddLog(const char* msg, const char* lbl) {
    auto old_size = Buf.size();
    Buf.appendf("[%s] %s\n", lbl, msg);
    for (auto new_size = Buf.size(); old_size < new_size; old_size++)
      if (Buf[old_size] == '\n') LineOffsets.push_back(old_size);
    ScrollToBottom = true;
  }

  void Draw() {
    if (ImGui::Button("Clear")) Clear();
    ImGui::SameLine();
    bool copy = ImGui::Button("Copy");
    ImGui::SameLine();
    Filter.Draw("Filter", -100.0f);
    ImGui::Separator();
    ImGui::BeginChild("scrolling");
    ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(0, 1));
    if (copy) ImGui::LogToClipboard();

    if (Filter.IsActive()) {
      const char* buf_begin = Buf.begin();
      const char* line      = buf_begin;
      for (int line_no = 0; line != nullptr; line_no++) {
        const char* line_end = (line_no < LineOffsets.Size)
                                   ? buf_begin + LineOffsets[line_no]
                                   : nullptr;
        if (Filter.PassFilter(line, line_end))
          ImGui::TextUnformatted(line, line_end);
        line = line_end != nullptr && line_end[1] != 0 ? line_end + 1 : nullptr;
      }
    } else {
      ImGui::TextUnformatted(Buf.begin());
    }

    if (ScrollToBottom) ImGui::SetScrollHere(1.0f);
    ScrollToBottom = false;
    ImGui::PopStyleVar();
    ImGui::EndChild();
  }
  void Draw(const char* title, bool* p_opened = nullptr) {
    ImGui::SetNextWindowSize(ImVec2(500, 400), ImGuiSetCond_FirstUseEver);
    ImGui::Begin(title, p_opened);
    Draw();
    ImGui::End();
  }
};

std::mutex  _log_mutex;
ImGuiAppLog _log_widget;
void        log_info(gui_window* win, const string& msg) {
  _log_mutex.lock();
  _log_widget.AddLog(msg.c_str(), "info");
  _log_mutex.unlock();
}
void log_error(gui_window* win, const string& msg) {
  _log_mutex.lock();
  _log_widget.AddLog(msg.c_str(), "errn");
  _log_mutex.unlock();
}
void clear_log(gui_window* win) {
  _log_mutex.lock();
  _log_widget.Clear();
  _log_mutex.unlock();
}
void draw_log(gui_window* win) {
  _log_mutex.lock();
  _log_widget.Draw();
  _log_mutex.unlock();
}

// draw param
bool draw_param(gui_window* win, const string& name, gui_param& param) {
  auto copy = param;
  switch (param.type) {
    case gui_param_type::value1f:
      if (param.minmaxf.x == param.minmaxf.y) {
        return draw_dragger(win, name.c_str(),
                   param.readonly ? (float&)copy.valuef
                                  : (float&)param.valuef) &&
               !param.readonly;
      } else {
        return draw_slider(win, name.c_str(),
                   param.readonly ? (float&)copy.valuef : (float&)param.valuef,
                   param.minmaxf.x, param.minmaxf.y) &&
               !param.readonly;
      }
      break;
    case gui_param_type::value2f:
      if (param.minmaxf.x == param.minmaxf.y) {
        return draw_dragger(win, name.c_str(),
                   param.readonly ? (vec2f&)copy.valuef
                                  : (vec2f&)param.valuef) &&
               !param.readonly;
      } else {
        return draw_slider(win, name.c_str(),
                   param.readonly ? (vec2f&)copy.valuef : (vec2f&)param.valuef,
                   param.minmaxf.x, param.minmaxf.y) &&
               !param.readonly;
      }
      break;
    case gui_param_type::value3f:
      if (param.color) {
        return draw_coloredit(win, name.c_str(),
                   param.readonly ? (vec3f&)copy.valuef
                                  : (vec3f&)param.valuef) &&
               !param.readonly;
      } else if (param.minmaxf.x == param.minmaxf.y) {
        return draw_dragger(win, name.c_str(),
                   param.readonly ? (vec3f&)copy.valuef
                                  : (vec3f&)param.valuef) &&
               !param.readonly;
      } else {
        return draw_slider(win, name.c_str(),
                   param.readonly ? copy.valuef : param.valuef, param.minmaxf.x,
                   param.minmaxf.y) &&
               !param.readonly;
      }
      break;
    case gui_param_type::value4f:
      if (param.color) {
        return draw_coloredit(win, name.c_str(),
                   param.readonly ? (vec4f&)copy.valuef
                                  : (vec4f&)param.valuef) &&
               !param.readonly;
      } else if (param.minmaxf.x == param.minmaxf.y) {
        return draw_dragger(win, name.c_str(),
                   param.readonly ? (vec4f&)copy.valuef
                                  : (vec4f&)param.valuef) &&
               !param.readonly;
      } else {
        return draw_slider(win, name.c_str(),
                   param.readonly ? (vec4f&)copy.valuef : (vec4f&)param.valuef,
                   param.minmaxf.x, param.minmaxf.y) &&
               !param.readonly;
      }
      break;
    case gui_param_type::value1i:
      if (!param.labels.empty()) {
        return draw_combobox(win, name.c_str(),
                   param.readonly ? (int&)copy.valuei : (int&)param.valuei,
                   param.labels) &&
               !param.readonly;
      } else if (param.minmaxi.x == param.minmaxi.y) {
        return draw_dragger(win, name.c_str(),
                   param.readonly ? (int&)copy.valuei : (int&)param.valuei) &&
               !param.readonly;
      } else {
        return draw_slider(win, name.c_str(),
                   param.readonly ? (int&)copy.valuei : (int&)param.valuei,
                   param.minmaxi.x, param.minmaxi.y) &&
               !param.readonly;
      }
      break;
    case gui_param_type::value2i:
      if (param.minmaxi.x == param.minmaxi.y) {
        return draw_dragger(win, name.c_str(),
                   param.readonly ? (vec2i&)copy.valuei
                                  : (vec2i&)param.valuei) &&
               !param.readonly;
      } else {
        return draw_slider(win, name.c_str(),
                   param.readonly ? (vec2i&)copy.valuei : (vec2i&)param.valuei,
                   param.minmaxi.x, param.minmaxi.y) &&
               !param.readonly;
      }
      break;
    case gui_param_type::value3i:
      if (param.minmaxi.x == param.minmaxi.y) {
        return draw_dragger(win, name.c_str(),
                   param.readonly ? (vec3i&)copy.valuei
                                  : (vec3i&)param.valuei) &&
               !param.readonly;
      } else {
        return draw_slider(win, name.c_str(),
                   param.readonly ? (vec3i&)copy.valuei : (vec3i&)param.valuei,
                   param.minmaxi.x, param.minmaxi.y) &&
               !param.readonly;
      }
      break;
    case gui_param_type::value4i:
      if (param.minmaxi.x == param.minmaxi.y) {
        return draw_dragger(win, name.c_str(),
                   param.readonly ? (vec4i&)copy.valuei
                                  : (vec4i&)param.valuei) &&
               !param.readonly;
      } else {
        return draw_slider(win, name.c_str(),
                   param.readonly ? (vec4i&)copy.valuei : (vec4i&)param.valuei,
                   param.minmaxi.x, param.minmaxi.y) &&
               !param.readonly;
      }
      break;
    case gui_param_type::value1s:
      if (!param.labels.empty()) {
        return draw_combobox(win, name.c_str(),
                   param.readonly ? copy.values : param.values, param.labels) &&
               !param.readonly;
      } else {
        return draw_textinput(win, name.c_str(),
                   param.readonly ? copy.values : param.values) &&
               !param.readonly;
      }
      break;
    case gui_param_type::value1b:
      if (!param.labels.empty()) {
        // maybe we should implement something different here
        return draw_checkbox(win, name.c_str(),
                   param.readonly ? copy.valueb : param.valueb) &&
               !param.readonly;
      } else {
        return draw_checkbox(win, name.c_str(),
                   param.readonly ? copy.valueb : param.valueb) &&
               !param.readonly;
      }
      break;
  }
}

// draw params
bool draw_params(gui_window* win, const string& name, gui_params& params) {
  auto edited = false;
  if (begin_header(win, name.c_str())) {
    for (auto& [name, param] : params) {
      auto pedited = draw_param(win, name, param);
      edited       = edited || pedited;
    }
    end_header(win);
  }
  return edited;
}

}  // namespace yocto
