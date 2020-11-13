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

#include "yocto_imgui.h"

#include <yocto/yocto_commonio.h>

#include <mutex>
#include <unordered_map>

#include "ext/glad/glad.h"

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif
#include <GLFW/glfw3.h>

#include "ext/imgui/imgui.h"
#include "ext/imgui/imgui_impl_glfw.h"
#include "ext/imgui/imgui_impl_opengl3.h"
#include "ext/imgui/imgui_internal.h"
#include "yocto_window.h"

#ifdef _WIN32
#undef near
#undef far
#endif

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::mutex;
using std::unordered_map;

}  // namespace yocto

// -----------------------------------------------------------------------------
// OPENGL WIDGETS
// -----------------------------------------------------------------------------
namespace yocto {

gui_widgets create_imgui(gui_window* window) {
  // widgets
  ImGui::CreateContext();
  ImGui::GetIO().IniFilename       = nullptr;
  ImGui::GetStyle().WindowRounding = 0;
  ImGui_ImplGlfw_InitForOpenGL(window->win, true);
#ifndef __APPLE__
  ImGui_ImplOpenGL3_Init();
#else
  ImGui_ImplOpenGL3_Init("#version 330");
#endif
  ImGui::StyleColorsDark();

  auto widgets   = gui_widgets{};
  widgets.window = window;
  return widgets;
}

bool begin_imgui(gui_widgets* widgets, const string& name,
    const vec2i& position, const vec2i& size) {
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGlfw_NewFrame();
  ImGui::NewFrame();
  //  if (win->widgets_left) {
  ImGui::SetNextWindowPos({0, 0});
  ImGui::SetNextWindowSize({(float)size.x, (float)size.y});
  //  } else {
  //    ImGui::SetNextWindowPos({(float)(window.x - win->widgets_width), 0});
  //    ImGui::SetNextWindowSize({(float)win->widgets_width, (float)window.y});
  //  }
  ImGui::SetNextWindowCollapsed(false);
  ImGui::SetNextWindowBgAlpha(1);

  return (ImGui::Begin(name.c_str(), nullptr,
      ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
          ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse |
          ImGuiWindowFlags_NoSavedSettings));

  // if (ImGui::Begin(win->title.c_str(), nullptr,
  //         ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
  //             ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse |
  //             ImGuiWindowFlags_NoSavedSettings)) {
  //   win->widgets_cb(win, win->input);
  // }
  // ImGui::End();
  // ImGui::Render();
  // ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

void end_imgui(gui_widgets* widgets) {
  ImGui::End();
  ImGui::Render();
  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

bool is_active(gui_widgets* widgets) {
  auto io = &ImGui::GetIO();
  return io->WantTextInput || io->WantCaptureMouse || io->WantCaptureKeyboard;
}

bool begin_header(gui_widgets* widgets, const char* lbl) {
  if (!ImGui::CollapsingHeader(lbl)) return false;
  ImGui::PushID(lbl);
  return true;
}
void end_header(gui_widgets* widgets) { ImGui::PopID(); }

void open_glmodal(gui_widgets* widgets, const char* lbl) {
  ImGui::OpenPopup(lbl);
}
void clear_glmodal(gui_widgets* widgets) { ImGui::CloseCurrentPopup(); }
bool begin_glmodal(gui_widgets* widgets, const char* lbl) {
  return ImGui::BeginPopupModal(lbl);
}
void end_glmodal(gui_widgets* widgets) { ImGui::EndPopup(); }
bool is_glmodal_open(gui_widgets* widgets, const char* lbl) {
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

bool draw_filedialog(gui_widgets* widgets, const char* lbl, string& path,
    bool save, const string& dirname, const string& filename,
    const string& filter) {
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

bool draw_filedialog_button(gui_widgets* widgets, const char* button_lbl,
    bool button_active, const char* lbl, string& path, bool save,
    const string& dirname, const string& filename, const string& filter) {
  if (is_glmodal_open(widgets, lbl)) {
    draw_button(widgets, button_lbl, button_active);
    return draw_filedialog(widgets, lbl, path, save, dirname, filename, filter);
  } else {
    if (draw_button(widgets, button_lbl, button_active)) {
      open_glmodal(widgets, lbl);
    }
    return false;
  }
}

bool draw_button(gui_widgets* widgets, const char* lbl, bool enabled) {
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

void draw_label(gui_widgets* widgets, const char* lbl, const string& label) {
  ImGui::LabelText(lbl, "%s", label.c_str());
}

void draw_separator(gui_widgets* widgets) { ImGui::Separator(); }

void continue_line(gui_widgets* widgets) { ImGui::SameLine(); }

bool draw_textinput(gui_widgets* widgets, const char* lbl, string& value) {
  auto buffer = array<char, 4096>{};
  auto num    = 0;
  for (auto c : value) buffer[num++] = c;
  buffer[num] = 0;
  auto edited = ImGui::InputText(lbl, buffer.data(), buffer.size());
  if (edited) value = buffer.data();
  return edited;
}

bool draw_slider(
    gui_widgets* widgets, const char* lbl, float& value, float min, float max) {
  return ImGui::SliderFloat(lbl, &value, min, max);
}
bool draw_slider(
    gui_widgets* widgets, const char* lbl, vec2f& value, float min, float max) {
  return ImGui::SliderFloat2(lbl, &value.x, min, max);
}
bool draw_slider(
    gui_widgets* widgets, const char* lbl, vec3f& value, float min, float max) {
  return ImGui::SliderFloat3(lbl, &value.x, min, max);
}
bool draw_slider(
    gui_widgets* widgets, const char* lbl, vec4f& value, float min, float max) {
  return ImGui::SliderFloat4(lbl, &value.x, min, max);
}

bool draw_slider(
    gui_widgets* widgets, const char* lbl, int& value, int min, int max) {
  return ImGui::SliderInt(lbl, &value, min, max);
}
bool draw_slider(
    gui_widgets* widgets, const char* lbl, vec2i& value, int min, int max) {
  return ImGui::SliderInt2(lbl, &value.x, min, max);
}
bool draw_slider(
    gui_widgets* widgets, const char* lbl, vec3i& value, int min, int max) {
  return ImGui::SliderInt3(lbl, &value.x, min, max);
}
bool draw_slider(
    gui_widgets* widgets, const char* lbl, vec4i& value, int min, int max) {
  return ImGui::SliderInt4(lbl, &value.x, min, max);
}

bool draw_dragger(gui_widgets* widgets, const char* lbl, float& value,
    float speed, float min, float max) {
  return ImGui::DragFloat(lbl, &value, speed, min, max);
}
bool draw_dragger(gui_widgets* widgets, const char* lbl, vec2f& value,
    float speed, float min, float max) {
  return ImGui::DragFloat2(lbl, &value.x, speed, min, max);
}
bool draw_dragger(gui_widgets* widgets, const char* lbl, vec3f& value,
    float speed, float min, float max) {
  return ImGui::DragFloat3(lbl, &value.x, speed, min, max);
}
bool draw_dragger(gui_widgets* widgets, const char* lbl, vec4f& value,
    float speed, float min, float max) {
  return ImGui::DragFloat4(lbl, &value.x, speed, min, max);
}

bool draw_dragger(gui_widgets* widgets, const char* lbl, int& value,
    float speed, int min, int max) {
  return ImGui::DragInt(lbl, &value, speed, min, max);
}
bool draw_dragger(gui_widgets* widgets, const char* lbl, vec2i& value,
    float speed, int min, int max) {
  return ImGui::DragInt2(lbl, &value.x, speed, min, max);
}
bool draw_dragger(gui_widgets* widgets, const char* lbl, vec3i& value,
    float speed, int min, int max) {
  return ImGui::DragInt3(lbl, &value.x, speed, min, max);
}
bool draw_dragger(gui_widgets* widgets, const char* lbl, vec4i& value,
    float speed, int min, int max) {
  return ImGui::DragInt4(lbl, &value.x, speed, min, max);
}

bool draw_checkbox(gui_widgets* widgets, const char* lbl, bool& value) {
  return ImGui::Checkbox(lbl, &value);
}
bool draw_checkbox(
    gui_widgets* widgets, const char* lbl, bool& value, bool invert) {
  if (!invert) {
    return draw_checkbox(widgets, lbl, value);
  } else {
    auto inverted = !value;
    auto edited   = ImGui::Checkbox(lbl, &inverted);
    if (edited) value = !inverted;
    return edited;
  }
}

bool draw_coloredit(gui_widgets* widgets, const char* lbl, vec3f& value) {
  auto flags = ImGuiColorEditFlags_Float;
  return ImGui::ColorEdit3(lbl, &value.x, flags);
}

bool draw_coloredit(gui_widgets* widgets, const char* lbl, vec4f& value) {
  auto flags = ImGuiColorEditFlags_Float;
  return ImGui::ColorEdit4(lbl, &value.x, flags);
}

bool draw_hdrcoloredit(gui_widgets* widgets, const char* lbl, vec3f& value) {
  auto color    = value;
  auto exposure = 0.0f;
  auto scale    = max(color);
  if (scale > 1) {
    color /= scale;
    exposure = log2(scale);
  }
  auto edit_exposure = draw_slider(
      widgets, (lbl + " [exp]"s).c_str(), exposure, 0, 10);
  auto edit_color = draw_coloredit(widgets, (lbl + " [col]"s).c_str(), color);
  if (edit_exposure || edit_color) {
    value = color * exp2(exposure);
    return true;
  } else {
    return false;
  }
}
bool draw_hdrcoloredit(gui_widgets* widgets, const char* lbl, vec4f& value) {
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
      widgets, (lbl + " [exp]"s).c_str(), exposure, 0, 10);
  auto edit_color = draw_coloredit(widgets, (lbl + " [col]"s).c_str(), color);
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

bool draw_combobox(gui_widgets* widgets, const char* lbl, int& value,
    const vector<string>& labels) {
  if (!ImGui::BeginCombo(lbl, labels[value].c_str())) return false;
  auto old_val = value;
  for (auto i = 0; i < labels.size(); i++) {
    ImGui::PushID(i);
    if (ImGui::Selectable(labels[i].c_str(), value == i)) value = i;
    if (value == i) ImGui::SetItemDefaultFocus();
    ImGui::PopID();
  }
  ImGui::EndCombo();
  return value != old_val;
}

bool draw_combobox(gui_widgets* widgets, const char* lbl, string& value,
    const vector<string>& labels) {
  if (!ImGui::BeginCombo(lbl, value.c_str())) return false;
  auto old_val = value;
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

bool draw_combobox(gui_widgets* widgets, const char* lbl, int& idx, int num,
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

void draw_progressbar(gui_widgets* widgets, const char* lbl, float fraction) {
  ImGui::PushStyleColor(ImGuiCol_PlotHistogram, ImVec4(0.5, 0.5, 1, 0.25));
  ImGui::ProgressBar(fraction, ImVec2(0.0f, 0.0f));
  ImGui::SameLine(0.0f, ImGui::GetStyle().ItemInnerSpacing.x);
  ImGui::Text(lbl, ImVec2(0.0f, 0.0f));
  ImGui::PopStyleColor(1);
}

void draw_progressbar(
    gui_widgets* widgets, const char* lbl, int current, int total) {
  auto overlay = std::to_string(current) + "/" + std::to_string(total);
  ImGui::PushStyleColor(ImGuiCol_PlotHistogram, ImVec4(0.5, 0.5, 1, 0.25));
  ImGui::ProgressBar(
      (float)current / (float)total, ImVec2(0.0f, 0.0f), overlay.c_str());
  ImGui::SameLine(0.0f, ImGui::GetStyle().ItemInnerSpacing.x);
  ImGui::Text(lbl, ImVec2(0.0f, 0.0f));
  ImGui::PopStyleColor(1);
}

void draw_histogram(
    gui_widgets* widgets, const char* lbl, const float* values, int count) {
  ImGui::PlotHistogram(lbl, values, count);
}
void draw_histogram(
    gui_widgets* widgets, const char* lbl, const vector<float>& values) {
  ImGui::PlotHistogram(lbl, values.data(), (int)values.size(), 0, nullptr,
      flt_max, flt_max, {0, 0}, 4);
}
void draw_histogram(
    gui_widgets* widgets, const char* lbl, const vector<vec2f>& values) {
  ImGui::PlotHistogram((lbl + " x"s).c_str(), (const float*)values.data() + 0,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec2f));
  ImGui::PlotHistogram((lbl + " y"s).c_str(), (const float*)values.data() + 1,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec2f));
}
void draw_histogram(
    gui_widgets* widgets, const char* lbl, const vector<vec3f>& values) {
  ImGui::PlotHistogram((lbl + " x"s).c_str(), (const float*)values.data() + 0,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec3f));
  ImGui::PlotHistogram((lbl + " y"s).c_str(), (const float*)values.data() + 1,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec3f));
  ImGui::PlotHistogram((lbl + " z"s).c_str(), (const float*)values.data() + 2,
      (int)values.size(), 0, nullptr, flt_max, flt_max, {0, 0}, sizeof(vec3f));
}
void draw_histogram(
    gui_widgets* widgets, const char* lbl, const vector<vec4f>& values) {
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
void        log_info(gui_widgets* widgets, const string& msg) {
  _log_mutex.lock();
  _log_widget.AddLog(msg.c_str(), "info");
  _log_mutex.unlock();
}
void log_error(gui_widgets* widgets, const string& msg) {
  _log_mutex.lock();
  _log_widget.AddLog(msg.c_str(), "errn");
  _log_mutex.unlock();
}
void clear_log(gui_widgets* widgets) {
  _log_mutex.lock();
  _log_widget.Clear();
  _log_mutex.unlock();
}
void draw_log(gui_widgets* widgets) {
  _log_mutex.lock();
  _log_widget.Draw();
  _log_mutex.unlock();
}

}  // namespace yocto
