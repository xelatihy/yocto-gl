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
#include "imgui/imgui_impl_glfw_gl3.h"

#include <map>

using uint = unsigned int;

inline uint make_glprogram(
    const char* vertex, const char* fragment, uint& vid, uint& fid, uint& vao) {
    assert(glGetError() == GL_NO_ERROR);
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    assert(glGetError() == GL_NO_ERROR);

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
    assert(glGetError() == GL_NO_ERROR);

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
    assert(glGetError() == GL_NO_ERROR);

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
    assert(glGetError() == GL_NO_ERROR);

    return pid;
}

inline uint make_glprogram(const char* vertex, const char* fragment) {
    uint vid = 0, fid = 0, vao = 0;
    return make_glprogram(vertex, fragment, vid, fid, vao);
}

// Create GLFW window
inline GLFWwindow* make_window(int width, int height, const std::string& title,
    void* user_ptr, void (*refresh)(GLFWwindow*)) {
    // glwindow
    if (!glfwInit()) throw std::runtime_error("cannot open glwindow");
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#if __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    auto win = glfwCreateWindow(width, height, title.c_str(), 0, 0);
    glfwMakeContextCurrent(win);
    glfwSwapInterval(1);  // Enable vsync

    glfwSetCharCallback(win, ImGui_ImplGlfw_CharCallback);
    glfwSetKeyCallback(win, ImGui_ImplGlfw_KeyCallback);
    glfwSetMouseButtonCallback(win, ImGui_ImplGlfw_MouseButtonCallback);
    glfwSetScrollCallback(win, ImGui_ImplGlfw_ScrollCallback);
    if (refresh) glfwSetWindowRefreshCallback(win, refresh);

    if (user_ptr) glfwSetWindowUserPointer(win, user_ptr);

// init gl extensions
#ifndef __APPLE__
    if (glewInit() != GLEW_OK) return nullptr;
#endif
    return win;
}

// Initialize ImGui widgets
inline void init_widgets(GLFWwindow* win) {
    ImGui::CreateContext();
    ImGui_ImplGlfwGL3_Init(win, false);
    ImGuiIO& io = ImGui::GetIO();
    // io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;  // Enable Keyboard
    // Controls io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;   // Enable
    io.IniFilename = nullptr;
}

// Begin draw widget
inline bool begin_widgets_frame(GLFWwindow* win, const char* title) {
    static auto first_time = true;
    ImGui_ImplGlfwGL3_NewFrame();
    if(first_time) {
        ImGui::SetNextWindowPos({0,0});
        ImGui::SetNextWindowSize({320,0});
        ImGui::SetNextWindowCollapsed(true);
        first_time = false;
    }
    return ImGui::Begin(title);
}

// End draw widget
inline void end_widgets_frame() {
    ImGui::End();
    ImGui::Render();
    ImGui_ImplGlfwGL3_RenderDrawData(ImGui::GetDrawData());
}

// Whether widget are active
inline bool get_widgets_active() {
    auto io = &ImGui::GetIO();
    return io->WantTextInput || io->WantCaptureMouse || io->WantCaptureKeyboard;
}

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
    for(auto c : *str) buf[num++] = c;
    buf[num] = 0;
    auto edited = InputText(label, buf, sizeof(buf));
    if(edited) *str = buf;
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
inline bool SelectableTreeNode(const char* lbl, std::shared_ptr<T>* selection,
    const std::shared_ptr<T>& content) {
    ImGuiTreeNodeFlags node_flags =
        ImGuiTreeNodeFlags_OpenOnArrow | ImGuiTreeNodeFlags_OpenOnDoubleClick;
    if (*selection == content) node_flags |= ImGuiTreeNodeFlags_Selected;
    auto open = ImGui::TreeNodeEx(content.get(), node_flags, "%s", lbl);
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
inline bool Combo(const char* lbl, std::shared_ptr<T>* val_,
    const std::vector<std::shared_ptr<T>>& vals, bool include_null) {
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
    scene_selection(const std::shared_ptr<T>& val)
        : ptr(val), tinfo(&typeid(T)) {}

    template <typename T>
    std::shared_ptr<T> as() const {
        return (&typeid(T) == tinfo) ? std::static_pointer_cast<T>(ptr) :
                                       nullptr;
    }

    std::shared_ptr<void> ptr = nullptr;    // selected pointer
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
    auto open = begin_glwidgets_tree(lbl.c_str(), selection, val);
    if (color != zero4f) ImGui::PopStyleColor();
    if (selection == val) sel = val;
    if (open) {
        draw_scene_tree_glwidgets_rec(lbl_, val, sel, highlights);
        ImGui::TreePop();
    }
}

template <typename T>
inline void draw_scene_tree_glwidgets_rec(const std::string& lbl_,
    const std::shared_ptr<T>& val, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {}

template <typename T>
inline void draw_glwidgets_scene_tree(const std::string& lbl_,
    const std::shared_ptr<T>& val, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    if (!val) return;
    auto lbl = val->name;
    if (!lbl_.empty()) lbl = lbl_ + ": " + val->name;
    auto selection = sel.as<T>();
    auto color = get_highlight_color(highlights, val->name);
    if (color != zero4f)
        ImGui::PushStyleColor(
            ImGuiCol_Text, {color.x, color.y, color.z, color.w});
    auto open = ImGui::SelectableTreeNode(lbl.c_str(), &selection, val);
    if (selection == val) sel = {val};
    if (color != zero4f) ImGui::PopStyleColor();
    if (selection == val) sel = val;
    if (open) {
        draw_scene_tree_glwidgets_rec(lbl_.c_str(), val, sel, highlights);
        ImGui::TreePop();
    }
}

template <>
inline void draw_scene_tree_glwidgets_rec<instance>(const std::string& lbl_,
    const std::shared_ptr<instance>& val, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    draw_glwidgets_scene_tree("shp", val->shp, sel, highlights);
    draw_glwidgets_scene_tree("sbd", val->sbd, sel, highlights);
    draw_glwidgets_scene_tree("mat", val->mat, sel, highlights);
}

template <>
inline void draw_scene_tree_glwidgets_rec<material>(const std::string& lbl_,
    const std::shared_ptr<material>& val, scene_selection& sel,
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
    const std::shared_ptr<environment>& val, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    draw_glwidgets_scene_tree("ke", val->ke_txt, sel, highlights);
}
template <>
inline void draw_scene_tree_glwidgets_rec<node>(const std::string& lbl_,
    const std::shared_ptr<node>& val, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    draw_glwidgets_scene_tree("ist", val->ist, sel, highlights);
    draw_glwidgets_scene_tree("cam", val->cam, sel, highlights);
    draw_glwidgets_scene_tree("env", val->env, sel, highlights);
    draw_glwidgets_scene_tree("par", val->parent, sel, highlights);
    auto cid = 0;
    for (auto ch : val->children) {
        draw_glwidgets_scene_tree(
            "ch" + std::to_string(cid++), ch.lock(), sel, highlights);
    }
}
template <>
inline void draw_scene_tree_glwidgets_rec<animation>(const std::string& lbl_,
    const std::shared_ptr<animation>& val, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    auto tid = 0;
    for (auto tg : val->targets) {
        draw_glwidgets_scene_tree(
            "tg" + std::to_string(tid++), tg, sel, highlights);
    }
}

inline void draw_glwidgets_scene_tree(const std::shared_ptr<scene>& scn,
    scene_selection& sel,
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
inline bool draw_glwidgets_scene_inspector(
    const std::shared_ptr<camera>& val, const std::shared_ptr<scene>& scn) {
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
inline bool draw_glwidgets_scene_inspector(
    const std::shared_ptr<texture>& val, const std::shared_ptr<scene>& scn) {
    auto edited = 0;
    edited += ImGui::InputText("name", &val->name);
    edited += ImGui::InputText("path", &val->path);
    edited += ImGui::Checkbox("clamp", &val->clamp);
    edited += ImGui::SliderFloat("scale", &val->scale, 0, 1);
    edited += ImGui::SliderFloat("gamma", &val->gamma, 1, 2.2f);
    ImGui::LabelText("img", "%d x %d", val->img.width(), val->img.height());
    return edited;
}

inline bool draw_glwidgets_scene_inspector(
    const std::shared_ptr<material>& val, const std::shared_ptr<scene>& scn) {
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

inline bool draw_glwidgets_scene_inspector(
    const std::shared_ptr<shape>& val, const std::shared_ptr<scene>& scn) {
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

inline bool draw_glwidgets_scene_inspector(
    const std::shared_ptr<subdiv>& val, const std::shared_ptr<scene>& scn) {
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

inline bool draw_glwidgets_scene_inspector(
    const std::shared_ptr<instance>& val, const std::shared_ptr<scene>& scn) {
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

inline bool draw_glwidgets_scene_inspector(
    const std::shared_ptr<environment>& val,
    const std::shared_ptr<scene>& scn) {
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

inline bool draw_glwidgets_scene_inspector(
    const std::shared_ptr<node>& val, const std::shared_ptr<scene>& scn) {
    auto edited = 0;
    edited += ImGui::InputText("name", &val->name);
    edited += ImGui::Combo("parent", &val->parent, scn->nodes, true);
    edited += ImGui::SliderFloat3("frame.x", &val->frame.x.x, -1, 1);
    edited += ImGui::SliderFloat3("frame.y", &val->frame.y.x, -1, 1);
    edited += ImGui::SliderFloat3("frame.z", &val->frame.z.x, -1, 1);
    edited += ImGui::SliderFloat3("frame.o", &val->frame.o.x, -10, 10);
    edited += ImGui::SliderFloat3("translation", &val->translation.x, -10, 10);
    edited += ImGui::SliderFloat4("rotation", &val->rotation.x, -1, 1);
    edited += ImGui::SliderFloat3("scale", &val->scale.x, 0, 10);
    edited += ImGui::Combo("cam", &val->cam, scn->cameras, true);
    edited += ImGui::Combo("ist", &val->ist, scn->instances, true);
    edited += ImGui::Combo("env", &val->env, scn->environments, true);
    return edited;
}

inline bool draw_glwidgets_scene_inspector(
    const std::shared_ptr<animation>& val, const std::shared_ptr<scene>& scn) {
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

inline bool draw_glwidgets_scene_tree(const std::string& lbl,
    const std::shared_ptr<scene>& scn, scene_selection& sel,
    std::vector<ygl::scene_selection>& update_list, int height,
    const std::unordered_map<std::string, std::string>& inspector_highlights) {
    if (!scn) return false;
    ImGui::PushID(scn.get());
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

inline bool draw_glwidgets_scene_inspector(const std::string& lbl,
    const std::shared_ptr<scene>& scn, scene_selection& sel,
    std::vector<ygl::scene_selection>& update_list, int height,
    const std::unordered_map<std::string, std::string>& inspector_highlights) {
    if (!scn || !sel.ptr) return false;
    ImGui::PushID(sel.ptr.get());
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
