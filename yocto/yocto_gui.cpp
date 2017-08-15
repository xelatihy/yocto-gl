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

// -----------------------------------------------------------------------------
// IMPLEMENTATION
// -----------------------------------------------------------------------------

#include "yocto_gui.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

// GL3W/GLFW
#ifdef __APPLE__
#define GLFW_INCLUDE_GLCOREARB
#include <GLFW/glfw3.h>
#include <OpenGL/gl3.h>
#else
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#endif

#include "ext/imgui/imgui.h"
#include "ext/imgui/imgui_impl_glfw_gl3.h"

namespace ygui {

//
// Window
//
struct window {
    GLFWwindow* win = nullptr;
    void* user_pointer = nullptr;

    int widget_width = 320;
    bool widget_enabled = false;

    text_callback text_cb = nullptr;
    mouse_callback mouse_cb = nullptr;
    refresh_callback refresh_cb = nullptr;
};

//
// Support
//
void glfw_error_cb(int error, const char* description) {
    printf("GLFW error: %s\n", description);
}

//
// Support
//
void glfw_text_cb(GLFWwindow* gwin, unsigned key) {
    auto win = (window*)glfwGetWindowUserPointer(gwin);
    if (win->widget_enabled) { ImGui_ImplGlfwGL3_CharCallback(win->win, key); }
    if (win->text_cb) win->text_cb(win, key);
}

//
// Support
//
void glfw_key_cb(
    GLFWwindow* gwin, int key, int scancode, int action, int mods) {
    auto win = (window*)glfwGetWindowUserPointer(gwin);
    if (win->widget_enabled) {
        ImGui_ImplGlfwGL3_KeyCallback(win->win, key, scancode, action, mods);
    }
}

//
// Support
//
void glfw_mouse_cb(GLFWwindow* gwin, int button, int action, int mods) {
    auto win = (window*)glfwGetWindowUserPointer(gwin);
    if (win->widget_enabled) {
        ImGui_ImplGlfwGL3_MouseButtonCallback(win->win, button, action, mods);
    }
    if (win->mouse_cb) win->mouse_cb(win, button, action == GLFW_PRESS, mods);
}

//
// Support
//
void glfw_scroll_cb(GLFWwindow* gwin, double xoffset, double yoffset) {
    auto win = (window*)glfwGetWindowUserPointer(gwin);
    if (win->widget_enabled) {
        ImGui_ImplGlfwGL3_ScrollCallback(win->win, xoffset, yoffset);
    }
}

//
// Support
//
void glfw_refresh_cb(GLFWwindow* gwin) {
    auto win = (window*)glfwGetWindowUserPointer(gwin);
    if (win->refresh_cb) win->refresh_cb(win);
}

//
// initialize glfw
//
window* init_window(int width, int height, const std::string& title,
    void* user_pointer, bool width_exclude_widgets) {
    // window
    auto win = new window();
    win->user_pointer = user_pointer;
    if (width_exclude_widgets) width += win->widget_width;

    // window
    if (!glfwInit()) return nullptr;

    // profile creation
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
#if __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    win->win = glfwCreateWindow(width, height, title.c_str(), 0, 0);
    glfwMakeContextCurrent(win->win);
    glfwSetWindowUserPointer(win->win, win);

    glfwSetErrorCallback(glfw_error_cb);

    glfwSetCharCallback(win->win, glfw_text_cb);
    glfwSetKeyCallback(win->win, glfw_key_cb);
    glfwSetMouseButtonCallback(win->win, glfw_mouse_cb);
    glfwSetScrollCallback(win->win, glfw_scroll_cb);

    glfwSetWindowRefreshCallback(win->win, glfw_refresh_cb);

// init gl extensions
#ifndef __APPLE__
    if (!glewInit()) return nullptr;
#endif

    return win;
}

//
// initialize glfw
//
void set_callbacks(window* win, text_callback text_cb, mouse_callback mouse_cb,
    refresh_callback refresh_cb) {
    win->text_cb = text_cb;
    win->mouse_cb = mouse_cb;
    win->refresh_cb = refresh_cb;
    if (text_cb) glfwSetCharCallback(win->win, glfw_text_cb);
}

//
// Clear glfw
//
void clear_window(window* win) {
    glfwDestroyWindow(win->win);
    glfwTerminate();
}

//
// Gets the user poiner
//
void* get_user_pointer(window* win) { return win->user_pointer; }

//
// Set window title
//
void set_window_title(window* win, const std::string& title) {
    glfwSetWindowTitle(win->win, title.c_str());
}

//
// Wait events
//
void wait_events(window* win) { glfwWaitEvents(); }

//
// Poll events
//
void poll_events(window* win) { glfwPollEvents(); }

//
// Swap buffers
//
void swap_buffers(window* win) { glfwSwapBuffers(win->win); }

//
// Should close
//
bool should_close(window* win) { return glfwWindowShouldClose(win->win); }

//
// Mouse button
//
int get_mouse_button(window* win) {
    auto mouse1 =
        glfwGetMouseButton(win->win, GLFW_MOUSE_BUTTON_1) == GLFW_PRESS;
    auto mouse2 =
        glfwGetMouseButton(win->win, GLFW_MOUSE_BUTTON_2) == GLFW_PRESS;
    auto mouse3 =
        glfwGetMouseButton(win->win, GLFW_MOUSE_BUTTON_3) == GLFW_PRESS;
    if (mouse1) return 1;
    if (mouse2) return 2;
    if (mouse3) return 3;
#if 0
            if (action == GLFW_RELEASE) {
                vparams.mouse_button = 0;
            } else if (button == GLFW_MOUSE_BUTTON_1 && !mods) {
                vparams.mouse_button = 1;
            } else if (button == GLFW_MOUSE_BUTTON_1 && (mods & GLFW_MOD_CONTROL)) {
                vparams.mouse_button = 2;
            } else if (button == GLFW_MOUSE_BUTTON_1 && (mods & GLFW_MOD_SHIFT)) {
                vparams.mouse_button = 3;
            } else if (button == GLFW_MOUSE_BUTTON_2) {
                vparams.mouse_button = 2;
            } else {
                vparams.mouse_button = 0;
            }
#endif
    return 0;
}

//
// Mouse position
//
ym::vec2i get_mouse_pos(window* win) {
    double x, y;
    glfwGetCursorPos(win->win, &x, &y);
    return {(int)x, (int)y};
}

//
// Mouse position
//
ym::vec2f get_mouse_posf(window* win) {
    double x, y;
    glfwGetCursorPos(win->win, &x, &y);
    return {(float)x, (float)y};
}

//
// Window size
//
ym::vec2i get_window_size(window* win) {
    auto ret = ym::vec2i{0, 0};
    glfwGetWindowSize(win->win, &ret[0], &ret[1]);
    return ret;
}

//
// Check if a key is pressed (not all keys are supported)
//
bool get_key(window* win, int key) {
    key = std::toupper(key);
    return glfwGetKey(win->win, key) == GLFW_PRESS;
}

//
// Framebuffer size
//
ym::vec2i get_framebuffer_size(window* win) {
    auto ret = ym::vec2i{0, 0};
    glfwGetFramebufferSize(win->win, &ret[0], &ret[1]);
    return ret;
}

//
// Widgets
//
int get_widget_size(window* win) { return win->widget_width; }

//
// Read pixels
//
std::vector<ym::vec4b> get_screenshot(
    window* win, ym::vec2i& wh, bool flipy, bool back) {
    wh = get_framebuffer_size(win);
    auto pixels = std::vector<ym::vec4b>(wh[0] * wh[1]);
    glReadBuffer((back) ? GL_BACK : GL_FRONT);
    glReadPixels(0, 0, wh[0], wh[1], GL_RGBA, GL_UNSIGNED_BYTE, pixels.data());
    if (flipy) {
        std::vector<ym::vec4b> line(wh[0]);
        for (int j = 0; j < wh[1] / 2; j++) {
            memcpy(line.data(), pixels.data() + j * wh[0] * 4, wh[0] * 4);
            memcpy(pixels.data() + j * wh[0] * 4,
                pixels.data() + (wh[1] - 1 - j) * wh[0] * 4, wh[0] * 4);
            memcpy(pixels.data() + (wh[1] - 1 - j) * wh[0] * 4, line.data(),
                wh[0] * 4);
        }
    }
    return pixels;
}

//
// ImGui
//
void init_widgets(window* win) {
    ImGui_ImplGlfwGL3_Init(win->win, false);
    ImGui::GetStyle().WindowRounding = 0;
    ImGui::GetIO().IniFilename = nullptr;
    win->widget_enabled = true;
}

//
// Nuklear
//
void clear_widgets(window* win) {
    ImGui_ImplGlfwGL3_Shutdown();
    win->widget_enabled = false;
}

//
// Begin draw widget
//
bool begin_widgets(window* win, const std::string& title) {
    ImGui_ImplGlfwGL3_NewFrame();
    auto size = get_window_size(win);
    ImGui::SetNextWindowSize({(float)win->widget_width, (float)size[1]});
    ImGui::SetNextWindowPos({(float)(size[0] - win->widget_width), (float)0});
    auto flags = ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize;
    ImGui::Begin(title.c_str(), nullptr, flags);
    // ImGui::ShowTestWindow();
    // ImGui::ShowStyleEditor();
    return true;
}

//
// End draw widget
//
void end_widgets(window* win) {
    ImGui::End();
    ImGui::Render();
}

//
// Whether widget are active
//
bool get_widget_active(window* win) {
    if (!win->widget_enabled) return false;
    auto io = &ImGui::GetIO();
    return io->WantTextInput || io->WantCaptureMouse || io->WantCaptureKeyboard;
}

//
// Horizontal separator
//
void separator_widget(window* win) { ImGui::Separator(); }

//
// Indent widget
//
void indent_begin_widgets(window* win) { ImGui::Indent(); }

//
// Indent widget
//
void indent_end_widgets(window* win) { ImGui::Unindent(); }

//
// Continue line with next widget
//
void continue_line_widgets(window* win) { ImGui::SameLine(); }

//
// Label widget
//
void label_widget(window* win, const std::string& lbl, const std::string& msg) {
    ImGui::LabelText(lbl.c_str(), "%s", msg.c_str());
}

//
// Label widget
//
void label_widget(window* win, const std::string& lbl, const char*& msg) {
    ImGui::LabelText(lbl.c_str(), "%s", msg);
}

//
// Label and int widget
//
void label_widget(
    window* win, const std::string& lbl, int val, const char* fmt) {
    ImGui::LabelText(lbl.c_str(), fmt, val);
}

//
// Label and int widget
//
void label_widget(
    window* win, const std::string& lbl, ym::vec2i val, const char* fmt) {
    ImGui::LabelText(lbl.c_str(), fmt, val[0], val[1]);
}

//
// Label and int widget
//
void label_widget(
    window* win, const std::string& lbl, ym::vec3i val, const char* fmt) {
    ImGui::LabelText(lbl.c_str(), fmt, val[0], val[1], val[2]);
}

//
// Label and int widget
//
void label_widget(
    window* win, const std::string& lbl, ym::vec4i val, const char* fmt) {
    ImGui::LabelText(lbl.c_str(), fmt, val[0], val[1], val[2], val[3]);
}

//
// Label and float widget
//
void label_widget(
    window* win, const std::string& lbl, float val, const char* fmt) {
    ImGui::LabelText(lbl.c_str(), fmt, val);
}

//
// Label and float widget
//
void label_widget(
    window* win, const std::string& lbl, ym::vec2f val, const char* fmt) {
    ImGui::LabelText(lbl.c_str(), fmt, val[0], val[1]);
}

//
// Label and float widget
//
void label_widget(
    window* win, const std::string& lbl, ym::vec3f val, const char* fmt) {
    ImGui::LabelText(lbl.c_str(), fmt, val[0], val[1], val[2]);
}

//
// Label and float widget
//
void label_widget(
    window* win, const std::string& lbl, ym::vec4f val, const char* fmt) {
    ImGui::LabelText(lbl.c_str(), fmt, val[0], val[1], val[2], val[3]);
}

//
// Text widget
//
bool text_widget(window* win, const std::string& lbl, char* buf, int buf_size) {
    return ImGui::InputText(lbl.c_str(), buf, buf_size);
}

//
// Label widget
//
bool slider_widget(
    window* win, const std::string& lbl, int* val, int min, int max, int incr) {
    return ImGui::SliderInt(lbl.c_str(), val, min, max);
}

//
// Label widget
//
bool slider_widget(window* win, const std::string& lbl, ym::vec2i* val, int min,
    int max, int incr) {
    return ImGui::SliderInt2(lbl.c_str(), (int*)val, min, max);
}

//
// Label widget
//
bool slider_widget(window* win, const std::string& lbl, ym::vec3i* val, int min,
    int max, int incr) {
    return ImGui::SliderInt3(lbl.c_str(), (int*)val, min, max);
}

//
// Label widget
//
bool slider_widget(window* win, const std::string& lbl, ym::vec4i* val, int min,
    int max, int incr) {
    return ImGui::SliderInt4(lbl.c_str(), (int*)val, min, max);
}

//
// Label widget
//
bool slider_widget(window* win, const std::string& lbl, float* val, float min,
    float max, float incr) {
    return ImGui::SliderFloat(lbl.c_str(), (float*)val, min, max);
}

//
// Label widget
//
bool slider_widget(window* win, const std::string& lbl, ym::vec2f* val,
    float min, float max, float incr) {
    return ImGui::SliderFloat2(lbl.c_str(), (float*)val, min, max);
}

//
// Label widget
//
bool slider_widget(window* win, const std::string& lbl, ym::vec3f* val,
    float min, float max, float incr) {
    return ImGui::SliderFloat3(lbl.c_str(), (float*)val, min, max);
}

//
// Label widget
//
bool slider_widget(window* win, const std::string& lbl, ym::vec4f* val,
    float min, float max, float incr) {
    return ImGui::SliderFloat4(lbl.c_str(), (float*)val, min, max);
}

//
// Slider widget
//
bool slider_widget(window* win, const std::string& lbl, ym::mat4f* val,
    float min, float max, float incr) {
    auto mod = false;
    for (auto i = 0; i < 4; i++) {
        auto modr = slider_widget(win, lbl + "[" + std::to_string(i + 1) + "]",
            val->data() + i, min, max, incr);
        mod = mod || modr;
    }
    return mod;
}

//
// Slider widget
//
bool slider_widget(
    window* win, const std::string& lbl, ym::quat4f* val, float incr) {
    auto mod = slider_widget(win, lbl, (ym::vec4f*)val, -1, 1, incr);
    if (mod) *val = ym::normalize(*val);
    return mod;
}

//
// Color widget
//
bool color_widget(window* win, const std::string& lbl, ym::vec4f* val) {
    return ImGui::ColorEdit4(lbl.c_str(), (float*)val);
}

//
// Color widget
//
bool color_widget(window* win, const std::string& lbl, ym::vec4b* val) {
    auto valf = ImGui::ColorConvertU32ToFloat4(*(uint32_t*)val);
    if (ImGui::ColorEdit4(lbl.c_str(), &valf.x)) {
        auto valb = ImGui::ColorConvertFloat4ToU32(valf);
        *(uint32_t*)val = valb;
        return true;
    }
    return false;
}

//
// Enum Widget
//
bool __enum_widget_labels_int(void* data, int idx, const char** out) {
    auto labels = (std::vector<std::pair<std::string, int>>*)data;
    *out = labels->at(idx).first.c_str();
    return true;
}

//
// Enum widget
//
bool combo_widget(window* win, const std::string& lbl, int* val,
    const std::vector<std::pair<std::string, int>>& labels) {
    auto cur = -1;
    for (auto idx = 0; idx < labels.size(); idx++) {
        if (labels[idx].second == *val) cur = idx;
    }
    assert(cur >= 0);
    auto ok = ImGui::Combo(lbl.c_str(), &cur, __enum_widget_labels_int,
        (void*)&labels, (int)labels.size());
    *val = labels[cur].second;
    return ok;
}

//
// List widget
//
bool list_widget(window* win, const std::string& lbl, int* val,
    const std::vector<std::pair<std::string, int>>& labels) {
    auto cur = -1;
    for (auto idx = 0; idx < labels.size(); idx++) {
        if (labels[idx].second == *val) cur = idx;
    }
    assert(cur >= 0);
    auto ok = ImGui::ListBox(lbl.c_str(), &cur, __enum_widget_labels_int,
        (void*)&labels, (int)labels.size());
    *val = labels[cur].second;
    return ok;
}

//
// Enum Widget
//
bool __enum_widget_labels_ptr(void* data, int idx, const char** out) {
    auto labels = (std::vector<std::pair<std::string, int>>*)data;
    *out = labels->at(idx).first.c_str();
    return true;
}

//
// Enum widget
//
bool combo_widget(window* win, const std::string& lbl, void** val,
    const std::vector<std::pair<std::string, void*>>& labels) {
    auto cur = -1;
    for (auto idx = 0; idx < labels.size(); idx++) {
        if (labels[idx].second == *val) cur = idx;
    }
    assert(cur >= 0);
    auto ok = ImGui::Combo(lbl.c_str(), &cur, __enum_widget_labels_ptr,
        (void*)&labels, (int)labels.size());
    *val = labels[cur].second;
    return ok;
}

//
// Enum widget
//
bool list_widget(window* win, const std::string& lbl, void** val,
    const std::vector<std::pair<std::string, void*>>& labels) {
    auto cur = -1;
    for (auto idx = 0; idx < labels.size(); idx++) {
        if (labels[idx].second == *val) cur = idx;
    }
    assert(cur >= 0);
    auto ok = ImGui::ListBox(lbl.c_str(), &cur, __enum_widget_labels_ptr,
        (void*)&labels, (int)labels.size());
    *val = labels[cur].second;
    return ok;
}

//
// Bool widget
//
bool checkbox_widget(window* win, const std::string& lbl, bool* val) {
    return ImGui::Checkbox(lbl.c_str(), val);
}

//
// Button widget
//
bool button_widget(window* win, const std::string& lbl) {
    return ImGui::Button(lbl.c_str());
}

//
// Collapsible header
//
bool collapsing_header_widget(window* win, const std::string& lbl) {
    return ImGui::CollapsingHeader(lbl.c_str());
}

//
// Start tree node
//
bool tree_begin_widget(window* win, const std::string& lbl) {
    return ImGui::TreeNode(lbl.c_str());
}

//
// Collapsible header
//
void tree_end_widget(window* win) { ImGui::TreePop(); }

//
// Start selectable tree node
//
bool tree_begin_widget(
    window* win, const std::string& lbl, void** selection, void* content) {
    ImGuiTreeNodeFlags node_flags =
        ImGuiTreeNodeFlags_OpenOnArrow | ImGuiTreeNodeFlags_OpenOnDoubleClick;
    if (*selection == content) node_flags |= ImGuiTreeNodeFlags_Selected;
    auto open = ImGui::TreeNodeEx(content, node_flags, "%s", lbl.c_str());
    if (ImGui::IsItemClicked()) *selection = content;
    return open;
}

//
// Start selectable tree node
//
bool tree_begin_widget(window* win, const std::string& lbl, void** selection,
    void* content, const ym::vec4f& col) {
    ImGui::PushStyleColor(ImGuiCol_Text, {col.x, col.y, col.z, col.w});
    auto ret = tree_begin_widget(win, lbl, selection, content);
    ImGui::PopStyleColor();
    return ret;
}

//
// End selectable tree node
//
void tree_end_widget(window* win, void* content) { ImGui::TreePop(); }

//
// Selectable tree leaf node
//
void tree_leaf_widget(
    window* win, const std::string& lbl, void** selection, void* content) {
    ImGuiTreeNodeFlags node_flags =
        ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen;
    if (*selection == content) node_flags |= ImGuiTreeNodeFlags_Selected;
    ImGui::TreeNodeEx(content, node_flags, "%s", lbl.c_str());
    if (ImGui::IsItemClicked()) *selection = content;
}

//
// Selectable tree leaf node
//
void tree_leaf_widget(window* win, const std::string& lbl, void** selection,
    void* content, const ym::vec4f& col) {
    ImGui::PushStyleColor(ImGuiCol_Text, {col.x, col.y, col.z, col.w});
    tree_leaf_widget(win, lbl, selection, content);
    ImGui::PopStyleColor();
}

//
// Image widget
//
void image_widget(
    window* win, int tid, const ym::vec2i& size, const ym::vec2i& imsize) {
    auto w = ImGui::GetContentRegionAvailWidth();
    auto s = ym::vec2f{(float)size.x, (float)size.y};
    auto a = (float)imsize.x / (float)imsize.y;
    if (!s.x && !s.y) {
        s.x = w;
        s.y = w / a;
    } else if (s.x && !s.y) {
        s.y = s.x / a;
    } else if (!s.x && s.y) {
        s.x = s.y * a;
    } else {
        auto as = s.x / s.y;
        if (as / a > 1) {
            s.x = s.y * a;
        } else {
            s.y = s.x / a;
        }
    }
    if (s.x > w) {
        s.x = w;
        s.y = w / a;
    }
    ImGui::Image((void*)(size_t)tid, {s.x, s.y});
}

//
// Scroll region
//
void scroll_region_begin_widget(
    window* win, const std::string& lbl, int height, bool border) {
    ImGui::BeginChild(lbl.c_str(), ImVec2(0, height), border);
}

//
// Scroll region
//
void scroll_region_end_widget(window* win) { ImGui::EndChild(); }

//
// Scroll region
//
void scroll_region_here_widget(window* win) { ImGui::SetScrollHere(); }

//
// Group ids
//
void groupid_begin_widget(window* win, int gid) { ImGui::PushID(gid); }

//
// Group ids
//
void groupid_begin_widget(window* win, void* gid) { ImGui::PushID(gid); }

//
// Group ids
//
void groupid_end_widget(window* win) { ImGui::PopID(); }

//
// Text color
//
void text_color_begin_widget(window* win, const ym::vec4f& color) {
    ImGui::PushStyleColor(
        ImGuiCol_Text, {color[0], color[1], color[2], color[3]});
}

//
// Text color
//
void text_color_end_widget(window* win) { ImGui::PopStyleColor(); }

}  // namespace ygui
