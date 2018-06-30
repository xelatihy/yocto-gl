//
// # Yocto/GLUI: Tiny OpenGL utilities for writing simple interactive apps.
//
// Small set of utilities to support writing OpenGL 3.3, manage
// windows with GLFW and draw immediate-mode widgets with ImGui.
//
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

#ifndef _YGLUI_H_
#define _YGLUI_H_

#include "../yocto/ygl.h"

#ifdef __APPLE__
#include <OpenGL/gl3.h>
#define GLFW_INCLUDE_GLCOREARB
#else
#include <GL/glew.h>
#endif
#include <GLFW/glfw3.h>
#include "imgui/imgui.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"

#include <map>

// -----------------------------------------------------------------------------
// LIGHTWEIGHT OPENGL UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

using uint = unsigned int;

inline void check_glerror() { assert(glGetError() == GL_NO_ERROR); }

inline uint make_glprogram(
    const char* vertex, const char* fragment, uint& vid, uint& fid, uint& vao) {
    check_glerror();
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    check_glerror();

    int errflags;
    char errbuf[10000];

    // create vertex
    vid = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vid, 1, &vertex, NULL);
    glCompileShader(vid);
    glGetShaderiv(vid, GL_COMPILE_STATUS, &errflags);
    if (!errflags) {
        glGetShaderInfoLog(vid, 10000, 0, errbuf);
        throw std::runtime_error(
            std::string("shader not compiled\n\n") + errbuf);
    }
    check_glerror();

    // create fragment
    fid = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fid, 1, &fragment, NULL);
    glCompileShader(fid);
    glGetShaderiv(fid, GL_COMPILE_STATUS, &errflags);
    if (!errflags) {
        glGetShaderInfoLog(fid, 10000, 0, errbuf);
        throw std::runtime_error(
            std::string("shader not compiled\n\n") + errbuf);
    }
    check_glerror();

    // create program
    auto pid = glCreateProgram();
    glAttachShader(pid, vid);
    glAttachShader(pid, fid);
    glLinkProgram(pid);
    glValidateProgram(pid);
    glGetProgramiv(pid, GL_LINK_STATUS, &errflags);
    if (!errflags) {
        glGetProgramInfoLog(pid, 10000, 0, errbuf);
        throw std::runtime_error(
            std::string("program not linked\n\n") + errbuf);
    }
    glGetProgramiv(pid, GL_VALIDATE_STATUS, &errflags);
    if (!errflags) {
        glGetProgramInfoLog(pid, 10000, 0, errbuf);
        throw std::runtime_error(
            std::string("program not linked\n\n") + errbuf);
    }
    check_glerror();

    return pid;
}

inline uint make_glprogram(const char* vertex, const char* fragment) {
    uint vid = 0, fid = 0, vao = 0;
    return make_glprogram(vertex, fragment, vid, fid, vao);
}

template <typename T>
inline uint make_glbuffer(
    const std::vector<T>& data, bool elems, bool dynamic = false) {
    auto bid = (uint)0;
    auto target = (elems) ? GL_ELEMENT_ARRAY_BUFFER : GL_ARRAY_BUFFER;
    auto usage = (dynamic) ? GL_DYNAMIC_DRAW : GL_STATIC_DRAW;
    check_glerror();
    glGenBuffers(1, &bid);
    glBindBuffer(target, bid);
    glBufferData(target, sizeof(T) * data.size(), data.data(), usage);
    glBindBuffer(target, 0);
    check_glerror();
    return bid;
}

template <typename T>
inline void update_glbuffer(uint bid, const std::vector<T>& data, bool elems) {
    auto target = (elems) ? GL_ELEMENT_ARRAY_BUFFER : GL_ARRAY_BUFFER;
    check_glerror();
    glBindBuffer(target, bid);
    glBufferSubData(target, 0, sizeof(T) * data.size(), data.data());
    glBindBuffer(target, 0);
    check_glerror();
}

inline uint make_gltexture(const image4f& img, bool linear, bool mipmap) {
    auto tid = (uint)0;
    check_glerror();
    glGenTextures(1, &tid);
    glBindTexture(GL_TEXTURE_2D, tid);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img.width, img.height, 0, GL_RGBA,
        GL_FLOAT, img.pxl.data());
    if (!linear) {
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    } else if (!mipmap) {
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    } else {
        glGenerateMipmap(GL_TEXTURE_2D);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(
            GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    }
    glBindTexture(GL_TEXTURE_2D, 0);
    check_glerror();
    return tid;
}

inline void update_gltexture(
    uint tid, const image4f& img, int old_width, int old_height, bool mipmap) {
    check_glerror();
    glBindTexture(GL_TEXTURE_2D, tid);
    if (img.width != old_width || img.height != old_height) {
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img.width, img.height, 0,
            GL_RGBA, GL_FLOAT, img.pxl.data());
    } else {
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, img.width, img.height, GL_RGBA,
            GL_FLOAT, img.pxl.data());
    }
    if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, 0);
    check_glerror();
}

inline void set_gluniform(uint loc, int val) {
    check_glerror();
    glUniform1i(loc, val);
    check_glerror();
}
inline void set_gluniform(uint loc, float val) {
    check_glerror();
    glUniform1f(loc, val);
    check_glerror();
}
inline void set_gluniform(uint loc, const vec2f& val) {
    check_glerror();
    glUniform2f(loc, val.x, val.y);
    check_glerror();
}
inline void set_gluniform(uint loc, const vec3f& val) {
    check_glerror();
    glUniform3f(loc, val.x, val.y, val.z);
    check_glerror();
}
inline void set_gluniform(uint loc, const vec4f& val) {
    check_glerror();
    glUniform4f(loc, val.x, val.y, val.z, val.w);
    check_glerror();
}
inline void set_gluniform(uint loc, const mat4f& val) {
    check_glerror();
    glUniformMatrix4fv(loc, 1, false, &val.x.x);
    check_glerror();
}

template <typename T>
inline void set_gluniform(uint pid, const char* name, const T& val) {
    check_glerror();
    set_gluniform(glGetUniformLocation(pid, name), val);
    check_glerror();
}

inline void set_gluniform_texture(uint loc, uint tid, int unit) {
    check_glerror();
    glActiveTexture(GL_TEXTURE0 + unit);
    glBindTexture(GL_TEXTURE_2D, tid);
    check_glerror();
    glUniform1i(loc, unit);
    check_glerror();
}

inline void set_gluniform_texture(
    uint pid, const char* name, uint tid, int unit) {
    check_glerror();
    set_gluniform_texture(glGetUniformLocation(pid, name), tid, unit);
    check_glerror();
}

inline void set_glvertattrib(uint loc, uint bid, const vec2f& val) {
    check_glerror();
    if (bid) {
        glEnableVertexAttribArray(loc);
        glBindBuffer(GL_ARRAY_BUFFER, bid);
        glVertexAttribPointer(loc, 2, GL_FLOAT, false, 0, 0);
    } else {
        glVertexAttrib2f(loc, val.x, val.y);
    }
    check_glerror();
}

inline void set_glvertattrib(uint loc, uint bid, const vec3f& val) {
    check_glerror();
    if (bid) {
        glEnableVertexAttribArray(loc);
        glBindBuffer(GL_ARRAY_BUFFER, bid);
        glVertexAttribPointer(loc, 3, GL_FLOAT, false, 0, 0);
    } else {
        glVertexAttrib3f(loc, val.x, val.y, val.z);
    }
    check_glerror();
}

inline void set_glvertattrib(uint loc, uint bid, const vec4f& val) {
    check_glerror();
    if (bid) {
        glEnableVertexAttribArray(loc);
        glBindBuffer(GL_ARRAY_BUFFER, bid);
        glVertexAttribPointer(loc, 4, GL_FLOAT, false, 0, 0);
    } else {
        glVertexAttrib4f(loc, val.x, val.y, val.z, val.w);
    }
    check_glerror();
}

template <typename T>
inline void set_glvertattrib(
    uint pid, const char* name, uint bid, const T& val) {
    check_glerror();
    set_glvertattrib(glGetAttribLocation(pid, name), bid, val);
    check_glerror();
}

inline void draw_glpoints(uint bid, int num) {
    check_glerror();
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bid);
    glDrawElements(GL_POINTS, 1 * num, GL_UNSIGNED_INT, 0);
    check_glerror();
}
inline void draw_gllines(uint bid, int num) {
    check_glerror();
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bid);
    glDrawElements(GL_LINES, 2 * num, GL_UNSIGNED_INT, 0);
    check_glerror();
}
inline void draw_gltriangles(uint bid, int num) {
    check_glerror();
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bid);
    glDrawElements(GL_TRIANGLES, 3 * num, GL_UNSIGNED_INT, 0);
    check_glerror();
}

inline void bind_glprog(uint pid) {
    check_glerror();
    glUseProgram(pid);
    check_glerror();
}

static const char* draw_image_vertex =
    R"(
    #version 330

    layout(location = 0) in vec2 vert_pos;
    layout(location = 1) in vec2 vert_texcoord;

    out vec2 texcoord;

    void main() {
        texcoord = vert_texcoord.xy;
        vec2 upos = vert_pos;
        gl_Position = vec4(upos.x, upos.y, 0, 1);
    }

    )";

static const char* draw_image_fragment =
    R"(
    #version 330

    in vec2 texcoord;
    uniform sampler2D img;
    out vec4 color;

    void main() {
        color = texture(img,texcoord);
    }
    )";

inline void draw_glimage(
    const image4f& img, vec2f imcenter, float imscale, vec2i win_size) {
    static uint gl_prog = 0, gl_pbo = 0, gl_tbo = 0, gl_ebo = 0, gl_txt = 0;
    static int gl_width = 0, gl_height = 0;

    auto imsize = vec2i{img.width, img.height};
    auto pos_ = std::vector<vec2f>{{0, 0}, {0, 1}, {1, 1}, {1, 0}};
    auto texcoord = std::vector<vec2f>{{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    auto triangles = std::vector<vec3i>{{0, 1, 2}, {0, 2, 3}};

    if (!gl_prog)
        gl_prog = make_glprogram(draw_image_vertex, draw_image_fragment);
    if (!gl_pbo) gl_pbo = make_glbuffer(pos_, false);
    if (!gl_tbo) gl_tbo = make_glbuffer(texcoord, false);
    if (!gl_ebo) gl_ebo = make_glbuffer(triangles, true);

    if (!gl_txt) {
        gl_txt = make_gltexture(img, false, false);
    } else {
        update_gltexture(gl_txt, img, gl_width, gl_height, false);
        gl_width = img.width;
        gl_height = img.height;
    }

    bind_glprog(gl_prog);

    set_gluniform_texture(gl_prog, "img", gl_txt, 0);

    auto pos = std::vector<vec2f>{{-imsize.x / 2.0f, -imsize.y / 2.0f},
        {imsize.x / 2.0f, -imsize.y / 2.0f}, {imsize.x / 2.0f, imsize.y / 2.0f},
        {-imsize.x / 2.0f, imsize.y / 2.0f}};
    for (auto& p : pos) p = p * imscale + imcenter;
    for (auto& p : pos) {
        p.x = 2 * p.x / win_size.x - 1;
        p.y = 2 * p.y / win_size.y - 1;
        p.y = -p.y;
    }
    update_glbuffer(gl_pbo, pos, false);

    set_glvertattrib(gl_prog, "vert_pos", gl_pbo, zero2f);
    set_glvertattrib(gl_prog, "vert_texcoord", gl_tbo, zero2f);

    draw_gltriangles(gl_ebo, 4);

    glUseProgram(0);
}

inline void draw_glimage(GLFWwindow* win, const image4f& img, vec2f& imcenter,
    float& imscale, bool zoom_to_fit, const vec4f& background) {
    auto win_size = zero2i, framebuffer_size = zero2i;
    glfwGetWindowSize(win, &win_size.x, &win_size.y);
    glfwGetFramebufferSize(win, &framebuffer_size.x, &framebuffer_size.y);

    center_image(
        imcenter, imscale, {img.width, img.height}, win_size, zoom_to_fit);

    check_glerror();
    glViewport(0, 0, framebuffer_size.x, framebuffer_size.y);
    glClearColor(background.x, background.y, background.z, background.w);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    check_glerror();

    draw_glimage(img, imcenter, imscale, win_size);

    check_glerror();
    glDisable(GL_BLEND);
    check_glerror();
}

// Create GLFW window
inline GLFWwindow* make_window(int w, int h, const std::string& title,
    void* user_ptr, void (*refresh)(GLFWwindow*)) {
    // glwindow
    if (!glfwInit()) throw std::runtime_error("cannot open glwindow");
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#if __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    auto win = glfwCreateWindow(w, h, title.c_str(), nullptr, nullptr);
    glfwMakeContextCurrent(win);
    glfwSwapInterval(1);  // Enable vsync

    // init gl extensions
#ifndef __APPLE__
    if (glewInit() != GLEW_OK) return nullptr;
#endif

    //    glfwSetCharCallback(win, ImGui_ImplGlfw_CharCallback);
    //    glfwSetKeyCallback(win, ImGui_ImplGlfw_KeyCallback);
    //    glfwSetMouseButtonCallback(win, ImGui_ImplGlfw_MouseButtonCallback);
    //    glfwSetScrollCallback(win, ImGui_ImplGlfw_ScrollCallback);

    if (refresh) glfwSetWindowRefreshCallback(win, refresh);
    if (user_ptr) glfwSetWindowUserPointer(win, user_ptr);

    return win;
}

// Initialize ImGui widgets
inline void init_widgets(GLFWwindow* win) {
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    (void)io;
    // io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;  // Enable Keyboard
    // Controls io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;   // Enable
    // Gamepad Controls
    io.IniFilename = nullptr;

    // ImGui_ImplGlfwGL3_Init(win, false);
    ImGui_ImplGlfw_InitForOpenGL(win, true);
    ImGui_ImplOpenGL3_Init();
    ImGui::StyleColorsDark();
}

// Begin draw widget
inline bool begin_widgets_frame(
    GLFWwindow* win, const char* title, bool* open) {
    static auto first_time = true;
    // Start the ImGui frame
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    if (first_time) {
        ImGui::SetNextWindowPos({0, 0});
        ImGui::SetNextWindowSize({320, 0});
        ImGui::SetNextWindowCollapsed(true);
        first_time = false;
    }
    return ImGui::Begin(title, open);
}

// End draw widget
inline void end_widgets_frame() {
    ImGui::End();
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

// Whether widget are active
inline bool get_widgets_active() {
    auto io = &ImGui::GetIO();
    return io->WantTextInput || io->WantCaptureMouse || io->WantCaptureKeyboard;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// GLFW EXTENSIONS
// -----------------------------------------------------------------------------

// GLFW extension
inline void glfwGetCursorPosExt(GLFWwindow* win, float* x, float* y) {
    auto xx = 0.0, yy = 0.0;
    glfwGetCursorPos(win, &xx, &yy);
    *x = (float)xx;
    *y = (float)yy;
}

// GLFW extension
inline int glfwGetMouseButtonIndexExt(GLFWwindow* win) {
    if (glfwGetMouseButton(win, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) return 1;
    if (glfwGetMouseButton(win, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS)
        return 2;
    if (glfwGetMouseButton(win, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_PRESS)
        return 3;
    return 0;
}

// GLFW extension
inline bool glfwGetAltKeyExt(GLFWwindow* win) {
    return glfwGetKey(win, GLFW_KEY_LEFT_ALT) == GLFW_PRESS ||
           glfwGetKey(win, GLFW_KEY_RIGHT_ALT) == GLFW_PRESS;
}
inline bool glfwGetCtrlKeyExt(GLFWwindow* win) {
    return glfwGetKey(win, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS ||
           glfwGetKey(win, GLFW_KEY_RIGHT_CONTROL) == GLFW_PRESS;
}
inline bool glfwGetShiftKeyExt(GLFWwindow* win) {
    return glfwGetKey(win, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
           glfwGetKey(win, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS;
}

// -----------------------------------------------------------------------------
// IMGUI EXTENSIONS
// -----------------------------------------------------------------------------

// ImGui extensions
namespace ImGui {
// Check if active widgets
inline bool GetWidgetsActiveExt() {
    auto io = &ImGui::GetIO();
    return io->WantTextInput || io->WantCaptureMouse || io->WantCaptureKeyboard;
}

// Input text
inline bool InputText(const char* label, std::string* str) {
    char buf[4096];
    auto num = 0;
    for (auto c : *str) buf[num++] = c;
    buf[num] = 0;
    auto edited = InputText(label, buf, sizeof(buf));
    if (edited) *str = buf;
    return edited;
}

// Start selectable tree node
inline bool SelectableTreeNode(
    const char* lbl, void** selection, void* content) {
    ImGuiTreeNodeFlags node_flags =
        ImGuiTreeNodeFlags_OpenOnArrow | ImGuiTreeNodeFlags_OpenOnDoubleClick;
    if (*selection == content) node_flags |= ImGuiTreeNodeFlags_Selected;
    auto open = ImGui::TreeNodeEx(content, node_flags, "%s", lbl);
    if (ImGui::IsItemClicked()) *selection = content;
    return open;
}

// Start selectable tree node
template <typename T>
inline bool SelectableTreeNode(
    const char* lbl, std::shared_ptr<T>* selection, T* content) {
    ImGuiTreeNodeFlags node_flags =
        ImGuiTreeNodeFlags_OpenOnArrow | ImGuiTreeNodeFlags_OpenOnDoubleClick;
    if (*selection == content) node_flags |= ImGuiTreeNodeFlags_Selected;
    auto open = ImGui::TreeNodeEx(content, node_flags, "%s", lbl);
    if (ImGui::IsItemClicked()) *selection = content;
    return open;
}

// Selectable tree leaf node
inline void SelectableTreeLeaf(
    const char* lbl, void*& selection, void* content) {
    ImGuiTreeNodeFlags node_flags =
        ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen;
    if (selection == content) node_flags |= ImGuiTreeNodeFlags_Selected;
    ImGui::TreeNodeEx(content, node_flags, "%s", lbl);
    if (ImGui::IsItemClicked()) selection = content;
}

// Combo widget.
inline bool Combo(const char* lbl, int* idx_, int nitems,
    const std::function<const char*(int)>& label) {
    auto& idx = *idx_;
    if (!ImGui::BeginCombo(lbl, label(idx))) return false;
    auto old_idx = idx;
    for (auto i = 0; i < nitems; i++) {
        ImGui::PushID(i);
        if (ImGui::Selectable(label(i), idx == i)) idx = i;
        if (idx == i) ImGui::SetItemDefaultFocus();
        ImGui::PopID();
    }
    ImGui::EndCombo();
    return idx != old_idx;
}

// Combo widget.
inline bool Combo(const char* lbl, std::string* val_,
    const std::vector<std::string>& labels) {
    auto& val = *val_;
    if (!ImGui::BeginCombo(lbl, val.c_str())) return false;
    auto old_val = val;
    for (auto i = 0; i < labels.size(); i++) {
        ImGui::PushID(i);
        if (ImGui::Selectable(labels[i].c_str(), val == labels[i]))
            val = labels[i];
        if (val == labels[i]) ImGui::SetItemDefaultFocus();
        ImGui::PopID();
    }
    ImGui::EndCombo();
    return val != old_val;
}

// Combo widget
template <typename T>
inline bool Combo(
    const char* lbl, T** val_, const std::vector<T*>& vals, bool include_null) {
    auto& val = *val_;
    if (!ImGui::BeginCombo(lbl, (val) ? val->name.c_str() : "<none>"))
        return false;
    auto old_val = val;
    if (include_null) {
        ImGui::PushID(100000);
        if (ImGui::Selectable("<none>", val == nullptr)) val = nullptr;
        if (val == nullptr) ImGui::SetItemDefaultFocus();
        ImGui::PopID();
    }
    for (auto i = 0; i < vals.size(); i++) {
        ImGui::PushID(i);
        if (ImGui::Selectable(vals[i]->name.c_str(), val == vals[i]))
            val = vals[i];
        if (val == vals[i]) ImGui::SetItemDefaultFocus();
        ImGui::PopID();
    }
    ImGui::EndCombo();
    return val != old_val;
}

}  // namespace ImGui

// Yocto/GL scene widgets
namespace ygl {

// Scene selection.
struct scene_selection {
    scene_selection() : ptr(nullptr), tinfo(nullptr) {}
    template <typename T>
    scene_selection(T* val) : ptr(val), tinfo(&typeid(T)) {}

    template <typename T>
    T* as() const {
        return (&typeid(T) == tinfo) ? (T*)ptr : nullptr;
    }

    void* ptr = nullptr;                    // selected pointer
    const std::type_info* tinfo = nullptr;  // type info
};

inline const std::map<animation_type, std::string>& animation_type_names() {
    static auto names = std::map<animation_type, std::string>{
        {animation_type::linear, "linear"},
        {animation_type::step, "step"},
        {animation_type::bezier, "bezier"},
    };
    return names;
}

inline vec4f get_highlight_color(
    const std::unordered_map<std::string, std::string>& highlights,
    const std::string& name) {
    static const std::unordered_map<std::string, vec4f>
        draw_visitor_highlight_colors = {{"red", {1, 0.5f, 0.5f, 1}},
            {"green", {0.5f, 1, 0.5f, 1}}, {"blue", {0.5f, 0.5f, 1, 1}}};

    if (!highlights.empty() && highlights.find(name) != highlights.end()) {
        return draw_visitor_highlight_colors.at(highlights.at(name));
    }
    return zero4f;
}

template <typename T>
inline void draw_scene_tree_glwidgets_rec(const std::string& lbl_, T* val,
    scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {}

template <typename T>
inline void draw_glwidgets_scene_tree(const std::string& lbl_, T* val,
    scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    if (!val) return;
    auto lbl = val->name;
    if (!lbl_.empty()) lbl = lbl_ + ": " + val->name;
    auto selection = sel.ptr;
    auto color = get_highlight_color(highlights, val->name);
    if (color != zero4f)
        ImGui::PushStyleColor(
            ImGuiCol_Text, {color.x, color.y, color.z, color.w});
    auto open = ImGui::SelectableTreeNode(lbl.c_str(), &selection, val);
    if (color != zero4f) ImGui::PopStyleColor();
    if (selection == val) sel = val;
    if (open) {
        draw_scene_tree_glwidgets_rec(lbl_, val, sel, highlights);
        ImGui::TreePop();
    }
}

template <>
inline void draw_scene_tree_glwidgets_rec<instance>(const std::string& lbl_,
    instance* val, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    draw_glwidgets_scene_tree("shp", val->shp, sel, highlights);
    draw_glwidgets_scene_tree("sbd", val->sbd, sel, highlights);
    draw_glwidgets_scene_tree("mat", val->mat, sel, highlights);
}

template <>
inline void draw_scene_tree_glwidgets_rec<material>(const std::string& lbl_,
    material* val, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    draw_glwidgets_scene_tree("ke", val->ke_txt, sel, highlights);
    draw_glwidgets_scene_tree("kd", val->kd_txt, sel, highlights);
    draw_glwidgets_scene_tree("ks", val->ks_txt, sel, highlights);
    draw_glwidgets_scene_tree("bump", val->bump_txt, sel, highlights);
    draw_glwidgets_scene_tree("disp", val->disp_txt, sel, highlights);
    draw_glwidgets_scene_tree("norm", val->norm_txt, sel, highlights);
}
template <>
inline void draw_scene_tree_glwidgets_rec<environment>(const std::string& lbl_,
    environment* val, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    draw_glwidgets_scene_tree("ke", val->ke_txt, sel, highlights);
}
template <>
inline void draw_scene_tree_glwidgets_rec<node>(const std::string& lbl_,
    node* val, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    draw_glwidgets_scene_tree("ist", val->ist, sel, highlights);
    draw_glwidgets_scene_tree("cam", val->cam, sel, highlights);
    draw_glwidgets_scene_tree("env", val->env, sel, highlights);
    draw_glwidgets_scene_tree("par", val->parent, sel, highlights);
    auto cid = 0;
    for (auto ch : val->children) {
        draw_glwidgets_scene_tree(
            "ch" + std::to_string(cid++), ch, sel, highlights);
    }
}
template <>
inline void draw_scene_tree_glwidgets_rec<animation>(const std::string& lbl_,
    animation* val, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    auto tid = 0;
    for (auto tg : val->targets) {
        draw_glwidgets_scene_tree(
            "tg" + std::to_string(tid++), tg, sel, highlights);
    }
}

inline void draw_glwidgets_scene_tree(scene* scn, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    if (!scn->cameras.empty() && ImGui::TreeNode("cameras")) {
        for (auto v : scn->cameras)
            draw_glwidgets_scene_tree("", v, sel, highlights);
        ImGui::TreePop();
    }
    if (!scn->shapes.empty() && ImGui::TreeNode("shapes")) {
        for (auto v : scn->shapes)
            draw_glwidgets_scene_tree("", v, sel, highlights);
        ImGui::TreePop();
    }
    if (!scn->subdivs.empty() && ImGui::TreeNode("subdivs")) {
        for (auto v : scn->subdivs)
            draw_glwidgets_scene_tree("", v, sel, highlights);
        ImGui::TreePop();
    }
    if (!scn->instances.empty() && ImGui::TreeNode("instances")) {
        for (auto v : scn->instances)
            draw_glwidgets_scene_tree("", v, sel, highlights);
        ImGui::TreePop();
    }
    if (!scn->materials.empty() && ImGui::TreeNode("materials")) {
        for (auto v : scn->materials)
            draw_glwidgets_scene_tree("", v, sel, highlights);
        ImGui::TreePop();
    }
    if (!scn->textures.empty() && ImGui::TreeNode("textures")) {
        for (auto v : scn->textures)
            draw_glwidgets_scene_tree("", v, sel, highlights);
        ImGui::TreePop();
    }
    if (!scn->environments.empty() && ImGui::TreeNode("environments")) {
        for (auto v : scn->environments)
            draw_glwidgets_scene_tree("", v, sel, highlights);
        ImGui::TreePop();
    }
    if (!scn->nodes.empty() && ImGui::TreeNode("nodes")) {
        for (auto v : scn->nodes)
            draw_glwidgets_scene_tree("", v, sel, highlights);
        ImGui::TreePop();
    }
    if (!scn->animations.empty() && ImGui::TreeNode("animations")) {
        for (auto v : scn->animations)
            draw_glwidgets_scene_tree("", v, sel, highlights);
        ImGui::TreePop();
    }
}

/// Visit struct elements.
inline bool draw_glwidgets_scene_inspector(camera* val, scene* scn) {
    auto edited = 0;
    edited += ImGui::InputText("name", &val->name);
    edited += ImGui::SliderFloat3("frame.x", &val->frame.x.x, -1, 1);
    edited += ImGui::SliderFloat3("frame.y", &val->frame.y.x, -1, 1);
    edited += ImGui::SliderFloat3("frame.z", &val->frame.z.x, -1, 1);
    edited += ImGui::SliderFloat3("frame.o", &val->frame.o.x, -10, 10);
    edited += ImGui::Checkbox("ortho", &val->ortho);
    edited += ImGui::SliderFloat("width", &val->width, 0.01f, 1);
    edited += ImGui::SliderFloat("height", &val->height, 0.01f, 1);
    edited += ImGui::SliderFloat("focal", &val->focal, 0.01f, 1);
    edited += ImGui::SliderFloat("focus", &val->focus, 0.01f, 1000);
    edited += ImGui::SliderFloat("aperture", &val->aperture, 0, 5);
    edited += ImGui::SliderFloat("near", &val->near, 0.01f, 10);
    edited += ImGui::SliderFloat("far", &val->far, 10, 10000);
    return edited;
}

/// Visit struct elements.
inline bool draw_glwidgets_scene_inspector(texture* val, scene* scn) {
    auto edited = 0;
    edited += ImGui::InputText("name", &val->name);
    edited += ImGui::InputText("path", &val->path);
    edited += ImGui::Checkbox("clamp", &val->clamp);
    edited += ImGui::SliderFloat("scale", &val->scale, 0, 1);
    edited += ImGui::SliderFloat("gamma", &val->gamma, 1, 2.2f);
    ImGui::LabelText("img", "%d x %d", val->img.width, val->img.height);
    return edited;
}

inline bool draw_glwidgets_scene_inspector(material* val, scene* scn) {
    auto edited = 0;
    edited += ImGui::InputText("name", &val->name);
    edited += ImGui::ColorEdit3("ke", &val->ke.x);  // TODO: HDR
    edited += ImGui::ColorEdit3("kd", &val->kd.x);
    edited += ImGui::ColorEdit3("ks", &val->ks.x);
    edited += ImGui::ColorEdit3("kt", &val->kt.x);
    edited += ImGui::SliderFloat("rs", &val->rs, 0, 1);
    edited += ImGui::SliderFloat("op", &val->op, 0, 1);
    edited += ImGui::Checkbox("fresnel", &val->fresnel);
    ImGui::SameLine();
    edited += ImGui::Checkbox("refract", &val->refract);
    edited += ImGui::Combo("ke txt", &val->ke_txt, scn->textures, true);
    edited += ImGui::Combo("kd txt", &val->kd_txt, scn->textures, true);
    edited += ImGui::Combo("ks txt", &val->ks_txt, scn->textures, true);
    edited += ImGui::Combo("kt txt", &val->kt_txt, scn->textures, true);
    edited += ImGui::Combo("op txt", &val->op_txt, scn->textures, true);
    edited += ImGui::Combo("rs txt", &val->rs_txt, scn->textures, true);
    edited += ImGui::Combo("bump txt", &val->bump_txt, scn->textures, true);
    edited += ImGui::Combo("disp txt", &val->disp_txt, scn->textures, true);
    edited += ImGui::Combo("norm txt", &val->norm_txt, scn->textures, true);
    edited += ImGui::Checkbox("base metallic", &val->base_metallic);
    edited += ImGui::Checkbox("glTF textures", &val->gltf_textures);
    return edited;
}

inline bool draw_glwidgets_scene_inspector(shape* val, scene* scn) {
    auto edited = 0;
    edited += ImGui::InputText("name", &val->name);
    edited += ImGui::InputText("path", &val->path);
    ImGui::LabelText("lines", "%ld", val->lines.size());
    ImGui::LabelText("triangles", "%ld", val->triangles.size());
    ImGui::LabelText("pos", "%ld", val->pos.size());
    ImGui::LabelText("norm", "%ld", val->norm.size());
    ImGui::LabelText("texcoord", "%ld", val->texcoord.size());
    ImGui::LabelText("color", "%ld", val->color.size());
    ImGui::LabelText("radius", "%ld", val->radius.size());
    ImGui::LabelText("tangsp", "%ld", val->tangsp.size());
    return edited;
}

inline bool draw_glwidgets_scene_inspector(subdiv* val, scene* scn) {
    auto edited = 0;
    edited += ImGui::InputText("name", &val->name);
    edited += ImGui::SliderInt("level", &val->level, 0, 10);
    edited += ImGui::Checkbox("catmull-clark", &val->catmull_clark);
    ImGui::SameLine();
    edited += ImGui::Checkbox("compute normals", &val->compute_normals);
    ImGui::LabelText("quads pos", "%ld", val->quads_pos.size());
    ImGui::LabelText("quads texcoord", "%ld", val->quads_texcoord.size());
    ImGui::LabelText("quads color", "%ld", val->quads_color.size());
    ImGui::LabelText("pos", "%ld", val->pos.size());
    ImGui::LabelText("texcoord", "%ld", val->texcoord.size());
    ImGui::LabelText("color", "%ld", val->color.size());
    return edited;
}

inline bool draw_glwidgets_scene_inspector(instance* val, scene* scn) {
    auto edited = 0;
    edited += ImGui::InputText("name", &val->name);
    edited += ImGui::SliderFloat3("frame.x", &val->frame.x.x, -1, 1);
    edited += ImGui::SliderFloat3("frame.y", &val->frame.y.x, -1, 1);
    edited += ImGui::SliderFloat3("frame.z", &val->frame.z.x, -1, 1);
    edited += ImGui::SliderFloat3("frame.o", &val->frame.o.x, -10, 10);
    edited += ImGui::Combo("shp", &val->shp, scn->shapes, true);
    edited += ImGui::Combo("sbd", &val->sbd, scn->subdivs, true);
    edited += ImGui::Combo("mat", &val->mat, scn->materials, true);
    return edited;
}

inline bool draw_glwidgets_scene_inspector(environment* val, scene* scn) {
    auto edited = 0;
    edited += ImGui::InputText("name", &val->name);
    edited += ImGui::SliderFloat3("frame.x", &val->frame.x.x, -1, 1);
    edited += ImGui::SliderFloat3("frame.y", &val->frame.y.x, -1, 1);
    edited += ImGui::SliderFloat3("frame.z", &val->frame.z.x, -1, 1);
    edited += ImGui::SliderFloat3("frame.o", &val->frame.o.x, -10, 10);
    edited += ImGui::ColorEdit4("ke", &val->ke.x);  // TODO: HDR
    edited += ImGui::Combo("ke txt", &val->ke_txt, scn->textures, true);
    return edited;
}

inline bool draw_glwidgets_scene_inspector(node* val, scene* scn) {
    auto edited = 0;
    edited += ImGui::InputText("name", &val->name);
    edited += ImGui::Combo("parent", &val->parent, scn->nodes, true);
    edited += ImGui::SliderFloat3("local.x", &val->local.x.x, -1, 1);
    edited += ImGui::SliderFloat3("local.y", &val->local.y.x, -1, 1);
    edited += ImGui::SliderFloat3("local.z", &val->local.z.x, -1, 1);
    edited += ImGui::SliderFloat3("local.o", &val->local.o.x, -10, 10);
    edited += ImGui::SliderFloat3("translation", &val->translation.x, -10, 10);
    edited += ImGui::SliderFloat4("rotation", &val->rotation.x, -1, 1);
    edited += ImGui::SliderFloat3("scale", &val->scale.x, 0, 10);
    edited += ImGui::Combo("cam", &val->cam, scn->cameras, true);
    edited += ImGui::Combo("ist", &val->ist, scn->instances, true);
    edited += ImGui::Combo("env", &val->env, scn->environments, true);
    return edited;
}

inline bool draw_glwidgets_scene_inspector(animation* val, scene* scn) {
    auto edited = 0;
    edited += ImGui::InputText("name", &val->name);
    edited += ImGui::InputText("path", &val->path);
    edited += ImGui::InputText("group", &val->group);
    // edited += ImGui::Combo("type", &val->type, animation_type_names());
    ImGui::LabelText("times", "%ld", val->times.size());
    ImGui::LabelText("translation", "%ld", val->translation.size());
    ImGui::LabelText("rotation", "%ld", val->rotation.size());
    ImGui::LabelText("scale", "%ld", val->scale.size());
    ImGui::LabelText("weights", "%ld", val->weights.size());
    ImGui::LabelText("targets", "%ld", val->targets.size());
    return edited;
}

inline bool draw_glwidgets_scene_tree(const std::string& lbl, scene* scn,
    scene_selection& sel, std::vector<ygl::scene_selection>& update_list,
    int height,
    const std::unordered_map<std::string, std::string>& inspector_highlights) {
    if (!scn) return false;
    ImGui::PushID(scn);
    ImGui::BeginChild("scrolling scene tree", ImVec2(0, height), false);
    draw_glwidgets_scene_tree(scn, sel, inspector_highlights);

    auto update_len = update_list.size();
#if 0
    if (test_scn) {
        draw_add_elem_glwidgets(
            scn, "cam", scn->cameras, test_scn->cameras, sel, update_list);
        draw_add_elem_glwidgets(scn, "txt", scn->textures,
            test_scn->textures, sel, update_list);
        draw_add_elem_glwidgets(scn, "mat", scn->materials,
            test_scn->materials, sel, update_list);
        draw_add_elem_glwidgets(
            scn, "shp", scn->shapes, test_scn->shapes, sel, update_list);
        draw_add_elem_glwidgets(scn, "ist", scn->instances,
            test_scn->instances, sel, update_list);
        draw_add_elem_glwidgets(
            scn, "nde", scn->nodes, test_scn->nodes, sel, update_list);
        draw_add_elem_glwidgets(scn, "env", scn->environments,
            test_scn->environments, sel, update_list);
        draw_add_elem_glwidgets(scn, "anim", scn->animations,
            test_scn->animations, sel, update_list);
    }
#endif

    ImGui::EndChild();
    ImGui::PopID();
    return update_list.size() != update_len;
}

inline bool draw_glwidgets_scene_inspector(const std::string& lbl, scene* scn,
    scene_selection& sel, std::vector<ygl::scene_selection>& update_list,
    int height,
    const std::unordered_map<std::string, std::string>& inspector_highlights) {
    if (!scn || !sel.ptr) return false;
    ImGui::PushID(sel.ptr);
    ImGui::BeginChild("scrolling scene inspector", ImVec2(0, height), false);

    auto update_len = update_list.size();

    auto edited = false;
    if (sel.as<camera>())
        edited = draw_glwidgets_scene_inspector(sel.as<camera>(), scn);
    if (sel.as<shape>())
        edited = draw_glwidgets_scene_inspector(sel.as<shape>(), scn);
    if (sel.as<subdiv>())
        edited = draw_glwidgets_scene_inspector(sel.as<subdiv>(), scn);
    if (sel.as<texture>())
        edited = draw_glwidgets_scene_inspector(sel.as<texture>(), scn);
    if (sel.as<material>())
        edited = draw_glwidgets_scene_inspector(sel.as<material>(), scn);
    if (sel.as<environment>())
        edited = draw_glwidgets_scene_inspector(sel.as<environment>(), scn);
    if (sel.as<instance>())
        edited = draw_glwidgets_scene_inspector(sel.as<instance>(), scn);
    if (sel.as<node>())
        edited = draw_glwidgets_scene_inspector(sel.as<node>(), scn);
    if (sel.as<animation>())
        edited = draw_glwidgets_scene_inspector(sel.as<animation>(), scn);
    if (edited) update_list.push_back(sel);

    ImGui::EndChild();
    ImGui::PopID();
    return update_list.size() != update_len;
}

}  // namespace ygl

#endif
