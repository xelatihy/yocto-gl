//
// Utilities to display a scene graph using ImGui.
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

#ifndef _YSCENEUI_H_
#define _YSCENEUI_H_

#include "../yocto/ygl.h"
#include "imgui/imgui.h"
#include "imgui/imgui_ext.h"
#ifdef _WIN32
#undef near
#undef far
#endif

#include <map>

inline const std::map<ygl::animation_type, std::string>& animation_type_names() {
    static auto names = std::map<ygl::animation_type, std::string>{
        {ygl::animation_type::linear, "linear"},
        {ygl::animation_type::step, "step"},
        {ygl::animation_type::bezier, "bezier"},
    };
    return names;
}

template <typename T>
inline void draw_scene_tree_glwidgets_rec(ygl::glwindow* win, const std::string& lbl_, T* val,
    void*& sel) {}

template <typename T>
inline void draw_glwidgets_scene_tree(ygl::glwindow* win, const std::string& lbl_, T* val,
    void*& sel) {
    if (!val) return;
    auto lbl = val->name;
    if (!lbl_.empty()) lbl = lbl_ + ": " + val->name;
    if(ygl::begin_imgui_selectabletreenode(win, lbl.c_str(), sel, val)) {
        draw_scene_tree_glwidgets_rec(win,lbl_, val, sel);
        ygl::end_imgui_treenode(win);
    }
}

template <>
inline void draw_scene_tree_glwidgets_rec<ygl::instance>(ygl::glwindow* win, const std::string& lbl_,
    ygl::instance* val, void*& sel) {
    draw_glwidgets_scene_tree(win,"shp", val->shp, sel);
    draw_glwidgets_scene_tree(win,"sbd", val->sbd, sel);
    draw_glwidgets_scene_tree(win,"mat", val->mat, sel);
}

template <>
inline void draw_scene_tree_glwidgets_rec<ygl::material>(ygl::glwindow* win, const std::string& lbl_,
    ygl::material* val, void*& sel) {
    draw_glwidgets_scene_tree(win,"ke", val->ke_txt, sel);
    draw_glwidgets_scene_tree(win,"kd", val->kd_txt, sel);
    draw_glwidgets_scene_tree(win,"ks", val->ks_txt, sel);
    draw_glwidgets_scene_tree(win,"bump", val->bump_txt, sel);
    draw_glwidgets_scene_tree(win,"disp", val->disp_txt, sel);
    draw_glwidgets_scene_tree(win,"norm", val->norm_txt, sel);
}
template <>
inline void draw_scene_tree_glwidgets_rec<ygl::environment>(ygl::glwindow* win, const std::string& lbl_,
    ygl::environment* val, void*& sel) {
    draw_glwidgets_scene_tree(win,"ke", val->ke_txt, sel);
}
template <>
inline void draw_scene_tree_glwidgets_rec<ygl::node>(ygl::glwindow* win, const std::string& lbl_,
    ygl::node* val, void*& sel) {
    draw_glwidgets_scene_tree(win,"ist", val->ist, sel);
    draw_glwidgets_scene_tree(win,"cam", val->cam, sel);
    draw_glwidgets_scene_tree(win,"env", val->env, sel);
    draw_glwidgets_scene_tree(win,"par", val->parent, sel);
    auto cid = 0;
    for (auto ch : val->children) {
        draw_glwidgets_scene_tree(win,
            "ch" + std::to_string(cid++), ch, sel);
    }
}
template <>
inline void draw_scene_tree_glwidgets_rec<ygl::animation>(ygl::glwindow* win, const std::string& lbl_,
    ygl::animation* val, void*& sel) {
    auto tid = 0;
    for (auto tg : val->targets) {
        draw_glwidgets_scene_tree(win,
            "tg" + std::to_string(tid++), tg, sel);
    }
}

inline void draw_glwidgets_scene_tree(ygl::glwindow* win, ygl::scene* scn, void*& sel) {
    if (!scn->cameras.empty() && ygl::begin_imgui_treenode(win, "cameras")) {
        for (auto v : scn->cameras)
            draw_glwidgets_scene_tree(win,"", v, sel);
        ygl::end_imgui_treenode(win);
    }
    if (!scn->shapes.empty() && ygl::begin_imgui_treenode(win, "shapes")) {
        for (auto v : scn->shapes)
            draw_glwidgets_scene_tree(win,"", v, sel);
        ygl::end_imgui_treenode(win);
    }
    if (!scn->subdivs.empty() && ygl::begin_imgui_treenode(win, "subdivs")) {
        for (auto v : scn->subdivs)
            draw_glwidgets_scene_tree(win,"", v, sel);
        ygl::end_imgui_treenode(win);
    }
    if (!scn->instances.empty() && ygl::begin_imgui_treenode(win, "instances")) {
        for (auto v : scn->instances)
            draw_glwidgets_scene_tree(win,"", v, sel);
        ygl::end_imgui_treenode(win);
    }
    if (!scn->materials.empty() && ygl::begin_imgui_treenode(win, "materials")) {
        for (auto v : scn->materials)
            draw_glwidgets_scene_tree(win,"", v, sel);
        ygl::end_imgui_treenode(win);
    }
    if (!scn->textures.empty() && ygl::begin_imgui_treenode(win, "textures")) {
        for (auto v : scn->textures)
            draw_glwidgets_scene_tree(win,"", v, sel);
        ygl::end_imgui_treenode(win);
    }
    if (!scn->environments.empty() && ygl::begin_imgui_treenode(win, "environments")) {
        for (auto v : scn->environments)
            draw_glwidgets_scene_tree(win,"", v, sel);
        ygl::end_imgui_treenode(win);
    }
    if (!scn->nodes.empty() && ygl::begin_imgui_treenode(win, "nodes")) {
        for (auto v : scn->nodes)
            draw_glwidgets_scene_tree(win,"", v, sel);
        ygl::end_imgui_treenode(win);
    }
    if (!scn->animations.empty() && ygl::begin_imgui_treenode(win, "animations")) {
        for (auto v : scn->animations)
            draw_glwidgets_scene_tree(win,"", v, sel);
        ygl::end_imgui_treenode(win);
    }
}

/// Visit struct elements.
inline bool draw_glwidgets_scene_inspector(ygl::glwindow* win, ygl::camera* val, ygl::scene* scn) {
    auto edited = 0;
    edited += ygl::draw_imgui_inputtext(win, "name", val->name);
    edited += ygl::draw_imgui_slider(win, "frame.x", val->frame.x.x, -1, 1);
    edited += ygl::draw_imgui_slider(win, "frame.y", val->frame.y.x, -1, 1);
    edited += ygl::draw_imgui_slider(win, "frame.z", val->frame.z.x, -1, 1);
    edited += ygl::draw_imgui_slider(win, "frame.o", val->frame.o.x, -10, 10);
    edited += ygl::draw_imgui_checkbox(win, "ortho", val->ortho);
    edited += ygl::draw_imgui_slider(win, "width", val->width, 0.01f, 1);
    edited += ygl::draw_imgui_slider(win, "height", val->height, 0.01f, 1);
    edited += ygl::draw_imgui_slider(win, "focal", val->focal, 0.01f, 1);
    edited += ygl::draw_imgui_slider(win, "focus", val->focus, 0.01f, 1000);
    edited += ygl::draw_imgui_slider(win, "aperture", val->aperture, 0, 5);
    edited += ygl::draw_imgui_slider(win, "near", val->near, 0.01f, 10);
    edited += ygl::draw_imgui_slider(win, "far", val->far, 10, 10000);
    return edited;
}

/// Visit struct elements.
inline bool draw_glwidgets_scene_inspector(ygl::glwindow* win, ygl::texture* val, ygl::scene* scn) {
    auto edited = 0;
    edited += ygl::draw_imgui_inputtext(win, "name", val->name);
    edited += ygl::draw_imgui_inputtext(win, "path", val->path);
    edited += ygl::draw_imgui_checkbox(win, "clamp", val->clamp);
    edited += ygl::draw_imgui_slider(win, "scale", val->scale, 0, 1);
    edited += ygl::draw_imgui_slider(win, "gamma", val->gamma, 1, 2.2f);
    ygl::draw_imgui_label(win, "img", "%d x %d", val->img.width, val->img.height);
    return edited;
}

inline bool draw_glwidgets_scene_inspector(ygl::glwindow* win, ygl::material* val, ygl::scene* scn) {
    auto edited = 0;
    edited += ygl::draw_imgui_inputtext(win, "name", val->name);
    edited += ygl::draw_imgui_coloredit(win, "ke", val->ke);  // TODO: HDR
    edited += ygl::draw_imgui_coloredit(win, "kd", val->kd);
    edited += ygl::draw_imgui_coloredit(win, "ks", val->ks);
    edited += ygl::draw_imgui_coloredit(win, "kt", val->kt);
    edited += ygl::draw_imgui_slider(win, "rs", val->rs, 0, 1);
    edited += ygl::draw_imgui_slider(win, "op", val->op, 0, 1);
    edited += ygl::draw_imgui_checkbox(win, "fresnel", val->fresnel);
    ygl::continue_imgui_line(win);
    edited += ygl::draw_imgui_checkbox(win, "refract", val->refract);
    edited += ygl::draw_imgui_coloredit(win, "vd", val->vd); // 0, 10
    edited += ygl::draw_imgui_coloredit(win, "va", val->va); // 0, 1
    edited += ygl::draw_imgui_slider(win, "vg", val->vg, -1, 1);
    edited += ygl::draw_imgui_combobox(win, "ke txt", val->ke_txt, scn->textures, true);
    edited += ygl::draw_imgui_combobox(win, "kd txt", val->kd_txt, scn->textures, true);
    edited += ygl::draw_imgui_combobox(win, "ks txt", val->ks_txt, scn->textures, true);
    edited += ygl::draw_imgui_combobox(win, "kt txt", val->kt_txt, scn->textures, true);
    edited += ygl::draw_imgui_combobox(win, "op txt", val->op_txt, scn->textures, true);
    edited += ygl::draw_imgui_combobox(win, "rs txt", val->rs_txt, scn->textures, true);
    edited += ygl::draw_imgui_combobox(win, "bump txt", val->bump_txt, scn->textures, true);
    edited += ygl::draw_imgui_combobox(win, "disp txt", val->disp_txt, scn->textures, true);
    edited += ygl::draw_imgui_combobox(win, "norm txt", val->norm_txt, scn->textures, true);
    edited += ygl::draw_imgui_checkbox(win, "base metallic", val->base_metallic);
    edited += ygl::draw_imgui_checkbox(win, "glTF textures", val->gltf_textures);
    return edited;
}

inline bool draw_glwidgets_scene_inspector(ygl::glwindow* win, ygl::shape* val, ygl::scene* scn) {
    auto edited = 0;
    edited += ygl::draw_imgui_inputtext(win, "name", val->name);
    edited += ygl::draw_imgui_inputtext(win, "path", val->path);
    ygl::draw_imgui_label(win, "lines", "%ld", val->lines.size());
    ygl::draw_imgui_label(win, "triangles", "%ld", val->triangles.size());
    ygl::draw_imgui_label(win, "pos", "%ld", val->pos.size());
    ygl::draw_imgui_label(win, "norm", "%ld", val->norm.size());
    ygl::draw_imgui_label(win, "texcoord", "%ld", val->texcoord.size());
    ygl::draw_imgui_label(win, "color", "%ld", val->color.size());
    ygl::draw_imgui_label(win, "radius", "%ld", val->radius.size());
    ygl::draw_imgui_label(win, "tangsp", "%ld", val->tangsp.size());
    return edited;
}

inline bool draw_glwidgets_scene_inspector(ygl::glwindow* win, ygl::subdiv* val, ygl::scene* scn) {
    auto edited = 0;
    edited += ygl::draw_imgui_inputtext(win, "name", val->name);
    edited += ygl::draw_imgui_slider(win, "level", val->level, 0, 10);
    edited += ygl::draw_imgui_checkbox(win, "catmull-clark", val->catmull_clark);
    ygl::continue_imgui_line(win);
    edited += ygl::draw_imgui_checkbox(win, "compute normals", val->compute_normals);
    ygl::draw_imgui_label(win, "quads pos", "%ld", val->quads_pos.size());
    ygl::draw_imgui_label(win, "quads texcoord", "%ld", val->quads_texcoord.size());
    ygl::draw_imgui_label(win, "quads color", "%ld", val->quads_color.size());
    ygl::draw_imgui_label(win, "pos", "%ld", val->pos.size());
    ygl::draw_imgui_label(win, "texcoord", "%ld", val->texcoord.size());
    ygl::draw_imgui_label(win, "color", "%ld", val->color.size());
    return edited;
}

inline bool draw_glwidgets_scene_inspector(ygl::glwindow* win, ygl::instance* val, ygl::scene* scn) {
    auto edited = 0;
    edited += ygl::draw_imgui_inputtext(win, "name", val->name);
    edited += ygl::draw_imgui_slider(win, "frame.x", val->frame.x, -1, 1);
    edited += ygl::draw_imgui_slider(win, "frame.y", val->frame.y, -1, 1);
    edited += ygl::draw_imgui_slider(win, "frame.z", val->frame.z, -1, 1);
    edited += ygl::draw_imgui_slider(win, "frame.o", val->frame.o, -10, 10);
    edited += ygl::draw_imgui_combobox(win, "shp", val->shp, scn->shapes, true);
    edited += ygl::draw_imgui_combobox(win, "sbd", val->sbd, scn->subdivs, true);
    edited += ygl::draw_imgui_combobox(win, "mat", val->mat, scn->materials, true);
    return edited;
}

inline bool draw_glwidgets_scene_inspector(ygl::glwindow* win, ygl::environment* val, ygl::scene* scn) {
    auto edited = 0;
    edited += ygl::draw_imgui_inputtext(win, "name", val->name);
    edited += ygl::draw_imgui_slider(win, "frame.x", val->frame.x, -1, 1);
    edited += ygl::draw_imgui_slider(win, "frame.y", val->frame.y, -1, 1);
    edited += ygl::draw_imgui_slider(win, "frame.z", val->frame.z, -1, 1);
    edited += ygl::draw_imgui_slider(win, "frame.o", val->frame.o, -10, 10);
    edited += ygl::draw_imgui_coloredit(win, "ke", val->ke);  // TODO: HDR
    edited += ygl::draw_imgui_combobox(win, "ke txt", val->ke_txt, scn->textures, true);
    return edited;
}

inline bool draw_glwidgets_scene_inspector(ygl::glwindow* win, ygl::node* val, ygl::scene* scn) {
    auto edited = 0;
    edited += ygl::draw_imgui_inputtext(win, "name", val->name);
    edited += ygl::draw_imgui_combobox(win, "parent", val->parent, scn->nodes, true);
    edited += ygl::draw_imgui_slider(win, "local.x", val->local.x, -1, 1);
    edited += ygl::draw_imgui_slider(win, "local.y", val->local.y, -1, 1);
    edited += ygl::draw_imgui_slider(win, "local.z", val->local.z, -1, 1);
    edited += ygl::draw_imgui_slider(win, "local.o", val->local.o, -10, 10);
    edited += ygl::draw_imgui_slider(win, "translation", val->translation, -10, 10);
    edited += ygl::draw_imgui_slider(win, "rotation", val->rotation, -1, 1);
    edited += ygl::draw_imgui_slider(win, "scale", val->scale, 0, 10);
    edited += ygl::draw_imgui_combobox(win, "cam", val->cam, scn->cameras, true);
    edited += ygl::draw_imgui_combobox(win, "ist", val->ist, scn->instances, true);
    edited += ygl::draw_imgui_combobox(win, "env", val->env, scn->environments, true);
    return edited;
}

inline bool draw_glwidgets_scene_inspector(ygl::glwindow* win, ygl::animation* val, ygl::scene* scn) {
    auto edited = 0;
    edited += ygl::draw_imgui_inputtext(win, "name", val->name);
    edited += ygl::draw_imgui_inputtext(win, "path", val->path);
    edited += ygl::draw_imgui_inputtext(win, "group", val->group);
    // edited += ygl::draw_imgui_combobox(win, "type", &val->type, animation_type_names());
    ygl::draw_imgui_label(win, "times", "%ld", val->times.size());
    ygl::draw_imgui_label(win, "translation", "%ld", val->translation.size());
    ygl::draw_imgui_label(win, "rotation", "%ld", val->rotation.size());
    ygl::draw_imgui_label(win, "scale", "%ld", val->scale.size());
    ygl::draw_imgui_label(win, "weights", "%ld", val->weights.size());
    ygl::draw_imgui_label(win, "targets", "%ld", val->targets.size());
    return edited;
}

inline bool draw_glwidgets_scene_tree(ygl::glwindow* win, const std::string& lbl, ygl::scene* scn,
    void*& sel, std::vector<std::pair<std::string, void*>>& update_list,
    int height) {
    if (!scn) return false;
    ImGui::PushID(scn);
    ImGui::BeginChild("scrolling scene tree", ImVec2(0, height), false);
    draw_glwidgets_scene_tree(win,scn, sel);

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

inline bool draw_glwidgets_scene_inspector(ygl::glwindow* win, const std::string& lbl, ygl::scene* scn,
    void*& sel, std::vector<std::pair<std::string, void*>>& update_list,
    int height) {
    if (!scn || !sel) return false;
    ImGui::PushID(sel);
    ImGui::BeginChild("scrolling scene inspector", ImVec2(0, height), false);

    auto update_len = update_list.size();

    for(auto cam : scn->cameras) {
        if(cam != sel) continue;
        if(draw_glwidgets_scene_inspector(win, cam, scn))
            update_list.push_back({"camera", cam});
    }
    for(auto shp : scn->shapes) {
        if(shp != sel) continue ;
        if(draw_glwidgets_scene_inspector(win, shp, scn))
            update_list.push_back({"shape", shp});
    }
    for(auto sbd : scn->subdivs) {
        if(sbd != sel) continue ;
        if(draw_glwidgets_scene_inspector(win, sbd, scn))
            update_list.push_back({"subdiv", sbd});
    }
    for(auto txt : scn->textures) {
        if(txt != sel) continue ;
        if(draw_glwidgets_scene_inspector(win, txt, scn))
            update_list.push_back({"texture", txt});
    }
    for(auto mat : scn->materials) {
        if(mat != sel) continue ;
        if(draw_glwidgets_scene_inspector(win, mat, scn))
            update_list.push_back({"material", mat});
    }
    for(auto env : scn->environments) {
        if(env != sel) continue ;
        if(draw_glwidgets_scene_inspector(win, env, scn))
            update_list.push_back({"environment", env});
    }
    for(auto ist : scn->instances) {
        if(ist != sel) continue ;
        if(draw_glwidgets_scene_inspector(win, ist, scn))
            update_list.push_back({"instance", ist});
    }
    for(auto nde : scn->nodes) {
        if(nde != sel) continue ;
        if(draw_glwidgets_scene_inspector(win, nde, scn))
            update_list.push_back({"node", nde});
    }
    for(auto anm : scn->animations) {
        if(anm != sel) continue ;
        if(draw_glwidgets_scene_inspector(win, anm, scn))
            update_list.push_back({"animation", anm});
    }

    ImGui::EndChild();
    ImGui::PopID();
    return update_list.size() != update_len;
}

#endif
