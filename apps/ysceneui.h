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
#include "yglutils.h"
using namespace ygl;

inline const map<animation_type, string>& animation_type_names() {
    static auto names = map<animation_type, string>{
        {animation_type::linear, "linear"},
        {animation_type::step, "step"},
        {animation_type::bezier, "bezier"},
    };
    return names;
}

template <typename T>
inline void draw_scene_tree_glwidgets_rec(
    glwindow* win, const string& lbl_, T* val, void*& sel) {}

template <typename T>
inline void draw_glwidgets_scene_tree(
    glwindow* win, const string& lbl_, T* val, void*& sel) {
    if (!val) return;
    auto lbl = val->name;
    if (!lbl_.empty()) lbl = lbl_ + ": " + val->name;
    if (begin_selectabletreenode_glwidget(win, lbl.c_str(), sel, val)) {
        draw_scene_tree_glwidgets_rec(win, lbl_, val, sel);
        end_treenode_glwidget(win);
    }
}

template <>
inline void draw_scene_tree_glwidgets_rec<instance>(
    glwindow* win, const string& lbl_, instance* val, void*& sel) {
    draw_glwidgets_scene_tree(win, "shp", val->shp, sel);
    draw_glwidgets_scene_tree(win, "sbd", val->sbd, sel);
    draw_glwidgets_scene_tree(win, "mat", val->mat, sel);
}

template <>
inline void draw_scene_tree_glwidgets_rec<material>(
    glwindow* win, const string& lbl_, material* val, void*& sel) {
    draw_glwidgets_scene_tree(win, "ke", val->ke_txt, sel);
    draw_glwidgets_scene_tree(win, "kd", val->kd_txt, sel);
    draw_glwidgets_scene_tree(win, "ks", val->ks_txt, sel);
    draw_glwidgets_scene_tree(win, "bump", val->bump_txt, sel);
    draw_glwidgets_scene_tree(win, "disp", val->disp_txt, sel);
    draw_glwidgets_scene_tree(win, "norm", val->norm_txt, sel);
}
template <>
inline void draw_scene_tree_glwidgets_rec<environment>(
    glwindow* win, const string& lbl_, environment* val, void*& sel) {
    draw_glwidgets_scene_tree(win, "ke", val->ke_txt, sel);
}
template <>
inline void draw_scene_tree_glwidgets_rec<node>(
    glwindow* win, const string& lbl_, node* val, void*& sel) {
    draw_glwidgets_scene_tree(win, "ist", val->ist, sel);
    draw_glwidgets_scene_tree(win, "cam", val->cam, sel);
    draw_glwidgets_scene_tree(win, "env", val->env, sel);
    draw_glwidgets_scene_tree(win, "par", val->parent, sel);
    auto cid = 0;
    for (auto ch : val->children) {
        draw_glwidgets_scene_tree(win, "ch" + to_string(cid++), ch, sel);
    }
}
template <>
inline void draw_scene_tree_glwidgets_rec<animation>(
    glwindow* win, const string& lbl_, animation* val, void*& sel) {
    auto tid = 0;
    for (auto tg : val->targets) {
        draw_glwidgets_scene_tree(win, "tg" + to_string(tid++), tg, sel);
    }
}

inline void draw_glwidgets_scene_tree(glwindow* win, scene* scn, void*& sel) {
    if (!scn->cameras.empty() && begin_treenode_glwidget(win, "cameras")) {
        for (auto v : scn->cameras) draw_glwidgets_scene_tree(win, "", v, sel);
        end_treenode_glwidget(win);
    }
    if (!scn->shapes.empty() && begin_treenode_glwidget(win, "shapes")) {
        for (auto v : scn->shapes) draw_glwidgets_scene_tree(win, "", v, sel);
        end_treenode_glwidget(win);
    }
    if (!scn->subdivs.empty() && begin_treenode_glwidget(win, "subdivs")) {
        for (auto v : scn->subdivs) draw_glwidgets_scene_tree(win, "", v, sel);
        end_treenode_glwidget(win);
    }
    if (!scn->instances.empty() && begin_treenode_glwidget(win, "instances")) {
        for (auto v : scn->instances)
            draw_glwidgets_scene_tree(win, "", v, sel);
        end_treenode_glwidget(win);
    }
    if (!scn->materials.empty() && begin_treenode_glwidget(win, "materials")) {
        for (auto v : scn->materials)
            draw_glwidgets_scene_tree(win, "", v, sel);
        end_treenode_glwidget(win);
    }
    if (!scn->textures.empty() && begin_treenode_glwidget(win, "textures")) {
        for (auto v : scn->textures) draw_glwidgets_scene_tree(win, "", v, sel);
        end_treenode_glwidget(win);
    }
    if (!scn->environments.empty() &&
        begin_treenode_glwidget(win, "environments")) {
        for (auto v : scn->environments)
            draw_glwidgets_scene_tree(win, "", v, sel);
        end_treenode_glwidget(win);
    }
    if (!scn->nodes.empty() && begin_treenode_glwidget(win, "nodes")) {
        for (auto v : scn->nodes) draw_glwidgets_scene_tree(win, "", v, sel);
        end_treenode_glwidget(win);
    }
    if (!scn->animations.empty() && begin_treenode_glwidget(win,
                                        "animation"
                                        "s")) {
        for (auto v : scn->animations)
            draw_glwidgets_scene_tree(win, "", v, sel);
        end_treenode_glwidget(win);
    }
}

/// Visit struct elements.
inline bool draw_glwidgets_scene_inspector(
    glwindow* win, camera* val, scene* scn) {
    auto edited = 0;
    edited += draw_textinput_glwidget(win, "name", val->name);
    edited += draw_slider_glwidget(win, "frame.x", val->frame.x.x, -1, 1);
    edited += draw_slider_glwidget(win, "frame.y", val->frame.y.x, -1, 1);
    edited += draw_slider_glwidget(win, "frame.z", val->frame.z.x, -1, 1);
    edited += draw_slider_glwidget(win, "frame.o", val->frame.o.x, -10, 10);
    edited += draw_checkbox_glwidget(win, "ortho", val->ortho);
    edited += draw_slider_glwidget(win, "film", val->film, 0.01f, 1);
    edited += draw_slider_glwidget(win, "focal", val->focal, 0.01f, 1);
    edited += draw_slider_glwidget(win, "focus", val->focus, 0.01f, 1000);
    edited += draw_slider_glwidget(win, "aperture", val->aperture, 0, 5);
    return edited;
}

/// Visit struct elements.
inline bool draw_glwidgets_scene_inspector(
    glwindow* win, texture* val, scene* scn) {
    auto edited = 0;
    edited += draw_textinput_glwidget(win, "name", val->name);
    edited += draw_textinput_glwidget(win, "path", val->path);
    edited += draw_checkbox_glwidget(win, "clamp", val->clamp);
    edited += draw_slider_glwidget(win, "scale", val->scale, 0, 1);
    edited += draw_checkbox_glwidget(win, "srgb", val->srgb);
    draw_label_glwidgets(
        win, "hdr", "%d x %d", width(val->imgf), height(val->imgf));
    draw_label_glwidgets(
        win, "ldr", "%d x %d", width(val->imgb), height(val->imgb));
    return edited;
}

inline bool draw_glwidgets_scene_inspector(
    glwindow* win, material* val, scene* scn) {
    auto edited = 0;
    edited += draw_textinput_glwidget(win, "name", val->name);
    edited += draw_coloredit_glwidget(win, "ke", val->ke);  // TODO: HDR
    edited += draw_coloredit_glwidget(win, "kd", val->kd);
    edited += draw_coloredit_glwidget(win, "ks", val->ks);
    edited += draw_coloredit_glwidget(win, "kt", val->kt);
    edited += draw_slider_glwidget(win, "rs", val->rs, 0, 1);
    edited += draw_slider_glwidget(win, "op", val->op, 0, 1);
    edited += draw_checkbox_glwidget(win, "fresnel", val->fresnel);
    continue_glwidgets_line(win);
    edited += draw_checkbox_glwidget(win, "refract", val->refract);
    edited += draw_coloredit_glwidget(win, "vd", val->vd);  // 0, 10
    edited += draw_coloredit_glwidget(win, "va", val->va);  // 0, 1
    edited += draw_slider_glwidget(win, "vg", val->vg, -1, 1);
    edited += draw_combobox_glwidget(
        win, "ke txt", val->ke_txt, scn->textures, true);
    edited += draw_combobox_glwidget(
        win, "kd txt", val->kd_txt, scn->textures, true);
    edited += draw_combobox_glwidget(
        win, "ks txt", val->ks_txt, scn->textures, true);
    edited += draw_combobox_glwidget(
        win, "kt txt", val->kt_txt, scn->textures, true);
    edited += draw_combobox_glwidget(
        win, "op txt", val->op_txt, scn->textures, true);
    edited += draw_combobox_glwidget(
        win, "rs txt", val->rs_txt, scn->textures, true);
    edited += draw_combobox_glwidget(
        win, "bump txt", val->bump_txt, scn->textures, true);
    edited += draw_combobox_glwidget(
        win, "disp txt", val->disp_txt, scn->textures, true);
    edited += draw_combobox_glwidget(
        win, "norm txt", val->norm_txt, scn->textures, true);
    edited += draw_checkbox_glwidget(win, "base metallic", val->base_metallic);
    edited += draw_checkbox_glwidget(win, "glTF textures", val->gltf_textures);
    return edited;
}

inline bool draw_glwidgets_scene_inspector(
    glwindow* win, shape* val, scene* scn) {
    auto edited = 0;
    edited += draw_textinput_glwidget(win, "name", val->name);
    edited += draw_textinput_glwidget(win, "path", val->path);
    draw_label_glwidgets(win, "lines", "%ld", val->lines.size());
    draw_label_glwidgets(win, "triangles", "%ld", val->triangles.size());
    draw_label_glwidgets(win, "pos", "%ld", val->pos.size());
    draw_label_glwidgets(win, "norm", "%ld", val->norm.size());
    draw_label_glwidgets(win, "texcoord", "%ld", val->texcoord.size());
    draw_label_glwidgets(win, "color", "%ld", val->color.size());
    draw_label_glwidgets(win, "radius", "%ld", val->radius.size());
    draw_label_glwidgets(win, "tangsp", "%ld", val->tangsp.size());
    return edited;
}

inline bool draw_glwidgets_scene_inspector(
    glwindow* win, subdiv* val, scene* scn) {
    auto edited = 0;
    edited += draw_textinput_glwidget(win, "name", val->name);
    edited += draw_slider_glwidget(win, "level", val->level, 0, 10);
    edited += draw_checkbox_glwidget(win, "catmull-clark", val->catmull_clark);
    continue_glwidgets_line(win);
    edited += draw_checkbox_glwidget(
        win, "compute normals", val->compute_normals);
    draw_label_glwidgets(win, "quads pos", "%ld", val->quads_pos.size());
    draw_label_glwidgets(
        win, "quads texcoord", "%ld", val->quads_texcoord.size());
    draw_label_glwidgets(win, "quads color", "%ld", val->quads_color.size());
    draw_label_glwidgets(win, "pos", "%ld", val->pos.size());
    draw_label_glwidgets(win, "texcoord", "%ld", val->texcoord.size());
    draw_label_glwidgets(win, "color", "%ld", val->color.size());
    return edited;
}

inline bool draw_glwidgets_scene_inspector(
    glwindow* win, instance* val, scene* scn) {
    auto edited = 0;
    edited += draw_textinput_glwidget(win, "name", val->name);
    edited += draw_slider_glwidget(win, "frame.x", val->frame.x, -1, 1);
    edited += draw_slider_glwidget(win, "frame.y", val->frame.y, -1, 1);
    edited += draw_slider_glwidget(win, "frame.z", val->frame.z, -1, 1);
    edited += draw_slider_glwidget(win, "frame.o", val->frame.o, -10, 10);
    edited += draw_combobox_glwidget(win, "shp", val->shp, scn->shapes, true);
    edited += draw_combobox_glwidget(win, "sbd", val->sbd, scn->subdivs, true);
    edited += draw_combobox_glwidget(win, "mat", val->mat, scn->materials, true);
    return edited;
}

inline bool draw_glwidgets_scene_inspector(
    glwindow* win, environment* val, scene* scn) {
    auto edited = 0;
    edited += draw_textinput_glwidget(win, "name", val->name);
    edited += draw_slider_glwidget(win, "frame.x", val->frame.x, -1, 1);
    edited += draw_slider_glwidget(win, "frame.y", val->frame.y, -1, 1);
    edited += draw_slider_glwidget(win, "frame.z", val->frame.z, -1, 1);
    edited += draw_slider_glwidget(win, "frame.o", val->frame.o, -10, 10);
    edited += draw_coloredit_glwidget(win, "ke", val->ke);  // TODO: HDR
    edited += draw_combobox_glwidget(
        win, "ke txt", val->ke_txt, scn->textures, true);
    return edited;
}

inline bool draw_glwidgets_scene_inspector(
    glwindow* win, node* val, scene* scn) {
    auto edited = 0;
    edited += draw_textinput_glwidget(win, "name", val->name);
    edited += draw_combobox_glwidget(
        win, "parent", val->parent, scn->nodes, true);
    edited += draw_slider_glwidget(win, "local.x", val->local.x, -1, 1);
    edited += draw_slider_glwidget(win, "local.y", val->local.y, -1, 1);
    edited += draw_slider_glwidget(win, "local.z", val->local.z, -1, 1);
    edited += draw_slider_glwidget(win, "local.o", val->local.o, -10, 10);
    edited += draw_slider_glwidget(
        win, "translation", val->translation, -10, 10);
    edited += draw_slider_glwidget(win, "rotation", val->rotation, -1, 1);
    edited += draw_slider_glwidget(win, "scale", val->scale, 0, 10);
    edited += draw_combobox_glwidget(win, "cam", val->cam, scn->cameras, true);
    edited += draw_combobox_glwidget(win, "ist", val->ist, scn->instances, true);
    edited += draw_combobox_glwidget(
        win, "env", val->env, scn->environments, true);
    return edited;
}

inline bool draw_glwidgets_scene_inspector(
    glwindow* win, animation* val, scene* scn) {
    auto edited = 0;
    edited += draw_textinput_glwidget(win, "name", val->name);
    edited += draw_textinput_glwidget(win, "path", val->path);
    edited += draw_textinput_glwidget(win, "group", val->group);
    // edited += draw_combobox_glwidget(win, "type", &val->type,
    // animation_type_names());
    draw_label_glwidgets(win, "times", "%ld", val->times.size());
    draw_label_glwidgets(win, "translation", "%ld", val->translation.size());
    draw_label_glwidgets(win, "rotation", "%ld", val->rotation.size());
    draw_label_glwidgets(win, "scale", "%ld", val->scale.size());
    draw_label_glwidgets(win, "weights", "%ld", val->weights.size());
    draw_label_glwidgets(win, "targets", "%ld", val->targets.size());
    return edited;
}

inline bool draw_glwidgets_scene_tree(glwindow* win, const string& lbl,
    scene* scn, void*& sel, vector<tuple<string, void*>>& update_list,
    int height) {
    if (!scn) return false;
    draw_glwidgets_scene_tree(win, scn, sel);
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
    return update_list.size() != update_len;
}

inline bool draw_glwidgets_scene_inspector(glwindow* win, const string& lbl,
    scene* scn, void*& sel, vector<tuple<string, void*>>& update_list,
    int height) {
    if (!scn || !sel) return false;
    begin_child_glwidget(win, "scrolling scene inspector", {0, height});

    auto update_len = update_list.size();

    for (auto cam : scn->cameras) {
        if (cam != sel) continue;
        if (draw_glwidgets_scene_inspector(win, cam, scn))
            update_list.push_back({"camera", cam});
    }
    for (auto shp : scn->shapes) {
        if (shp != sel) continue;
        if (draw_glwidgets_scene_inspector(win, shp, scn))
            update_list.push_back({"shape", shp});
    }
    for (auto sbd : scn->subdivs) {
        if (sbd != sel) continue;
        if (draw_glwidgets_scene_inspector(win, sbd, scn))
            update_list.push_back({"subdiv", sbd});
    }
    for (auto txt : scn->textures) {
        if (txt != sel) continue;
        if (draw_glwidgets_scene_inspector(win, txt, scn))
            update_list.push_back({"texture", txt});
    }
    for (auto mat : scn->materials) {
        if (mat != sel) continue;
        if (draw_glwidgets_scene_inspector(win, mat, scn))
            update_list.push_back({"material", mat});
    }
    for (auto env : scn->environments) {
        if (env != sel) continue;
        if (draw_glwidgets_scene_inspector(win, env, scn))
            update_list.push_back({"environment", env});
    }
    for (auto ist : scn->instances) {
        if (ist != sel) continue;
        if (draw_glwidgets_scene_inspector(win, ist, scn))
            update_list.push_back({"instance", ist});
    }
    for (auto nde : scn->nodes) {
        if (nde != sel) continue;
        if (draw_glwidgets_scene_inspector(win, nde, scn))
            update_list.push_back({"node", nde});
    }
    for (auto anm : scn->animations) {
        if (anm != sel) continue;
        if (draw_glwidgets_scene_inspector(win, anm, scn))
            update_list.push_back({"animation", anm});
    }

    end_child_glwidget(win);
    return update_list.size() != update_len;
}

#endif
