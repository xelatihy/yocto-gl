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

inline const map<yocto_interpolation_type, string>& animation_type_names() {
    static auto names = map<yocto_interpolation_type, string>{
        {yocto_interpolation_type::linear, "linear"},
        {yocto_interpolation_type::step, "step"},
        {yocto_interpolation_type::bezier, "bezier"},
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
inline void draw_scene_tree_glwidgets_rec<yocto_shape>(
    glwindow* win, const string& lbl_, yocto_shape* val, void*& sel) {
    draw_glwidgets_scene_tree(win, "material", val->material, sel);
}

template <>
inline void draw_scene_tree_glwidgets_rec<yocto_instance>(
    glwindow* win, const string& lbl_, yocto_instance* val, void*& sel) {
    draw_glwidgets_scene_tree(win, "shape", val->shape, sel);
}

template <>
inline void draw_scene_tree_glwidgets_rec<yocto_material>(
    glwindow* win, const string& lbl_, yocto_material* val, void*& sel) {
    draw_glwidgets_scene_tree(win, "emission", val->emission_texture, sel);
    draw_glwidgets_scene_tree(win, "diffuse", val->diffuse_texture, sel);
    draw_glwidgets_scene_tree(win, "specular", val->specular_texture, sel);
    draw_glwidgets_scene_tree(win, "bump", val->bump_texture, sel);
    draw_glwidgets_scene_tree(win, "displament", val->displacement_texture, sel);
    draw_glwidgets_scene_tree(win, "normal", val->normal_texture, sel);
}
template <>
inline void draw_scene_tree_glwidgets_rec<yocto_environment>(
    glwindow* win, const string& lbl_, yocto_environment* val, void*& sel) {
    draw_glwidgets_scene_tree(win, "emission", val->emission_texture, sel);
}
template <>
inline void draw_scene_tree_glwidgets_rec<yocto_scene_node>(
    glwindow* win, const string& lbl_, yocto_scene_node* val, void*& sel) {
    draw_glwidgets_scene_tree(win, "instance", val->instance, sel);
    draw_glwidgets_scene_tree(win, "camera", val->camera, sel);
    draw_glwidgets_scene_tree(win, "environment", val->environment, sel);
    draw_glwidgets_scene_tree(win, "parent", val->parent, sel);
    auto cid = 0;
    for (auto ch : val->children) {
        draw_glwidgets_scene_tree(win, "child" + to_string(cid++), ch, sel);
    }
}
template <>
inline void draw_scene_tree_glwidgets_rec<yocto_animation>(
    glwindow* win, const string& lbl_, yocto_animation* val, void*& sel) {
    auto tid = 0;
    for (auto tg : val->node_targets) {
        draw_glwidgets_scene_tree(win, "target" + to_string(tid++), tg, sel);
    }
}

inline void draw_glwidgets_scene_tree(
    glwindow* win, yocto_scene* scene, void*& sel) {
    if (!scene->cameras.empty() && begin_treenode_glwidget(win, "cameras")) {
        for (auto v : scene->cameras)
            draw_glwidgets_scene_tree(win, "", v, sel);
        end_treenode_glwidget(win);
    }
    if (!scene->shapes.empty() && begin_treenode_glwidget(win, "shapes")) {
        for (auto v : scene->shapes) draw_glwidgets_scene_tree(win, "", v, sel);
        end_treenode_glwidget(win);
    }
    if (!scene->instances.empty() && begin_treenode_glwidget(win, "instances")) {
        for (auto v : scene->instances)
            draw_glwidgets_scene_tree(win, "", v, sel);
        end_treenode_glwidget(win);
    }
    if (!scene->materials.empty() && begin_treenode_glwidget(win, "materials")) {
        for (auto v : scene->materials)
            draw_glwidgets_scene_tree(win, "", v, sel);
        end_treenode_glwidget(win);
    }
    if (!scene->textures.empty() && begin_treenode_glwidget(win, "textures")) {
        for (auto v : scene->textures)
            draw_glwidgets_scene_tree(win, "", v, sel);
        end_treenode_glwidget(win);
    }
    if (!scene->environments.empty() &&
        begin_treenode_glwidget(win, "environments")) {
        for (auto v : scene->environments)
            draw_glwidgets_scene_tree(win, "", v, sel);
        end_treenode_glwidget(win);
    }
    if (!scene->nodes.empty() && begin_treenode_glwidget(win, "nodes")) {
        for (auto v : scene->nodes) draw_glwidgets_scene_tree(win, "", v, sel);
        end_treenode_glwidget(win);
    }
    if (!scene->animations.empty() && begin_treenode_glwidget(win,
                                          "animation"
                                          "s")) {
        for (auto v : scene->animations)
            draw_glwidgets_scene_tree(win, "", v, sel);
        end_treenode_glwidget(win);
    }
}

/// Visit struct elements.
inline bool draw_glwidgets_scene_inspector(
    glwindow* win, yocto_camera* val, yocto_scene* scene) {
    auto edited = 0;
    edited += draw_textinput_glwidget(win, "name", val->name);
    edited += draw_slider_glwidget(win, "frame.x", val->frame.x.x, -1, 1);
    edited += draw_slider_glwidget(win, "frame.y", val->frame.y.x, -1, 1);
    edited += draw_slider_glwidget(win, "frame.z", val->frame.z.x, -1, 1);
    edited += draw_slider_glwidget(win, "frame.o", val->frame.o.x, -10, 10);
    edited += draw_checkbox_glwidget(win, "ortho", val->orthographic);
    edited += draw_slider_glwidget(win, "film", val->film_size, 0.01f, 1);
    edited += draw_slider_glwidget(win, "focal", val->focal_length, 0.01f, 1);
    edited += draw_slider_glwidget(
        win, "focus", val->focus_distance, 0.01f, 1000);
    edited += draw_slider_glwidget(win, "aperture", val->lens_aperture, 0, 5);
    return edited;
}

/// Visit struct elements.
inline bool draw_glwidgets_scene_inspector(
    glwindow* win, yocto_texture* val, yocto_scene* scene) {
    auto edited = 0;
    edited += draw_textinput_glwidget(win, "name", val->name);
    edited += draw_textinput_glwidget(win, "path", val->filename);
    edited += draw_checkbox_glwidget(win, "clamp_to_edge", val->clamp_to_edge);
    edited += draw_slider_glwidget(win, "scale", val->height_scale, 0, 1);
    edited += draw_checkbox_glwidget(win, "ldr_as_linear", val->ldr_as_linear);
    draw_label_glwidgets(win, "hdr_image", "%d x %d", val->hdr_image.width,
        val->hdr_image.height);
    draw_label_glwidgets(win, "ldr_image", "%d x %d", val->ldr_image.width,
        val->ldr_image.height);
    return edited;
}

inline bool draw_glwidgets_scene_inspector(
    glwindow* win, yocto_material* val, yocto_scene* scene) {
    auto edited = 0;
    edited += draw_textinput_glwidget(win, "name", val->name);
    edited += draw_coloredit_glwidget(win, "emission", val->emission);  // TODO:
                                                                        // HDR
    edited += draw_coloredit_glwidget(win, "diffuse", val->diffuse);
    edited += draw_coloredit_glwidget(win, "specular", val->specular);
    edited += draw_coloredit_glwidget(win, "transmission", val->transmission);
    edited += draw_slider_glwidget(win, "roughness", val->roughness, 0, 1);
    edited += draw_slider_glwidget(win, "opacity", val->opacity, 0, 1);
    edited += draw_checkbox_glwidget(win, "fresnel", val->fresnel);
    continue_glwidgets_line(win);
    edited += draw_checkbox_glwidget(win, "refract", val->refract);
    edited += draw_coloredit_glwidget(
        win, "volume_density", val->volume_density);  // 0, 10
    edited += draw_coloredit_glwidget(
        win, "volume_albedo", val->volume_albedo);  // 0, 1
    edited += draw_slider_glwidget(
        win, "volume_phaseg", val->volume_phaseg, -1, 1);
    edited += draw_combobox_glwidget(
        win, "emission_texture", val->emission_texture, scene->textures, true);
    edited += draw_combobox_glwidget(
        win, "diffuse_texture", val->diffuse_texture, scene->textures, true);
    edited += draw_combobox_glwidget(
        win, "specular_texture", val->specular_texture, scene->textures, true);
    edited += draw_combobox_glwidget(win, "transmission_texture",
        val->transmission_texture, scene->textures, true);
    edited += draw_combobox_glwidget(
        win, "opacity_texture", val->opacity_texture, scene->textures, true);
    edited += draw_combobox_glwidget(win, "roughness_texture",
        val->roughness_texture, scene->textures, true);
    edited += draw_combobox_glwidget(
        win, "bump_texture", val->bump_texture, scene->textures, true);
    edited += draw_combobox_glwidget(win, "displacement_texture",
        val->displacement_texture, scene->textures, true);
    edited += draw_combobox_glwidget(
        win, "normal_texture", val->normal_texture, scene->textures, true);
    edited += draw_checkbox_glwidget(win, "base metallic", val->base_metallic);
    edited += draw_checkbox_glwidget(win, "glTF textures", val->gltf_textures);
    return edited;
}

inline bool draw_glwidgets_scene_inspector(
    glwindow* win, yocto_shape* val, yocto_scene* scene) {
    auto edited = 0;
    edited += draw_textinput_glwidget(win, "name", val->name);
    edited += draw_textinput_glwidget(win, "path", val->filename);
    edited += draw_combobox_glwidget(
        win, "material", val->material, scene->materials, true);
    draw_label_glwidgets(win, "lines", "%ld", val->lines.size());
    draw_label_glwidgets(win, "triangles", "%ld", val->triangles.size());
    draw_label_glwidgets(win, "quads", "%ld", val->quads.size());
    draw_label_glwidgets(win, "pos", "%ld", val->positions.size());
    draw_label_glwidgets(win, "norm", "%ld", val->normals.size());
    draw_label_glwidgets(win, "texcoord", "%ld", val->texturecoords.size());
    draw_label_glwidgets(win, "color", "%ld", val->colors.size());
    draw_label_glwidgets(win, "radius", "%ld", val->radius.size());
    draw_label_glwidgets(win, "tangsp", "%ld", val->tangentspaces.size());
    return edited;
}

inline bool draw_glwidgets_scene_inspector(
    glwindow* win, yocto_instance* val, yocto_scene* scene) {
    auto edited = 0;
    edited += draw_textinput_glwidget(win, "name", val->name);
    edited += draw_slider_glwidget(win, "frame.x", val->frame.x, -1, 1);
    edited += draw_slider_glwidget(win, "frame.y", val->frame.y, -1, 1);
    edited += draw_slider_glwidget(win, "frame.z", val->frame.z, -1, 1);
    edited += draw_slider_glwidget(win, "frame.o", val->frame.o, -10, 10);
    edited += draw_combobox_glwidget(
        win, "shape", val->shape, scene->shapes, true);
    return edited;
}

inline bool draw_glwidgets_scene_inspector(
    glwindow* win, yocto_environment* val, yocto_scene* scene) {
    auto edited = 0;
    edited += draw_textinput_glwidget(win, "name", val->name);
    edited += draw_slider_glwidget(win, "frame.x", val->frame.x, -1, 1);
    edited += draw_slider_glwidget(win, "frame.y", val->frame.y, -1, 1);
    edited += draw_slider_glwidget(win, "frame.z", val->frame.z, -1, 1);
    edited += draw_slider_glwidget(win, "frame.o", val->frame.o, -10, 10);
    edited += draw_coloredit_glwidget(win, "ke", val->emission);  // TODO: HDR
    edited += draw_combobox_glwidget(
        win, "ke texture", val->emission_texture, scene->textures, true);
    return edited;
}

inline bool draw_glwidgets_scene_inspector(
    glwindow* win, yocto_scene_node* val, yocto_scene* scene) {
    auto edited = 0;
    edited += draw_textinput_glwidget(win, "name", val->name);
    edited += draw_combobox_glwidget(
        win, "parent", val->parent, scene->nodes, true);
    edited += draw_slider_glwidget(win, "local.x", val->local.x, -1, 1);
    edited += draw_slider_glwidget(win, "local.y", val->local.y, -1, 1);
    edited += draw_slider_glwidget(win, "local.z", val->local.z, -1, 1);
    edited += draw_slider_glwidget(win, "local.o", val->local.o, -10, 10);
    edited += draw_slider_glwidget(win, "translation", val->translation, -10, 10);
    edited += draw_slider_glwidget(win, "rotation", val->rotation, -1, 1);
    edited += draw_slider_glwidget(win, "scale", val->scale, 0, 10);
    edited += draw_combobox_glwidget(
        win, "camera", val->camera, scene->cameras, true);
    edited += draw_combobox_glwidget(
        win, "instance", val->instance, scene->instances, true);
    edited += draw_combobox_glwidget(
        win, "environment", val->environment, scene->environments, true);
    return edited;
}

inline bool draw_glwidgets_scene_inspector(
    glwindow* win, yocto_animation* val, yocto_scene* scene) {
    auto edited = 0;
    edited += draw_textinput_glwidget(win, "name", val->name);
    edited += draw_textinput_glwidget(win, "path", val->filename);
    edited += draw_textinput_glwidget(win, "group", val->animation_group);
    // edited += draw_combobox_glwidget(win, "type", &val->type,
    // animation_type_names());
    draw_label_glwidgets(win, "times", "%ld", val->keyframes_times.size());
    draw_label_glwidgets(
        win, "translation", "%ld", val->translation_keyframes.size());
    draw_label_glwidgets(win, "rotation", "%ld", val->rotation_keyframes.size());
    draw_label_glwidgets(win, "scale", "%ld", val->scale_keyframes.size());
    draw_label_glwidgets(
        win, "weights", "%ld", val->morph_weights_keyframes.size());
    draw_label_glwidgets(win, "targets", "%ld", val->node_targets.size());
    return edited;
}

inline bool draw_glwidgets_scene_tree(glwindow* win, const string& lbl,
    yocto_scene* scene, void*& sel, vector<tuple<string, void*>>& update_list,
    int height) {
    if (!scene) return false;
    draw_glwidgets_scene_tree(win, scene, sel);
    auto update_len = update_list.size();
#if 0
    if (test_scn) {
        draw_add_elem_glwidgets(
            scene, "camera", scene->cameras, test_scn->cameras, sel, update_list);
        draw_add_elem_glwidgets(scene, "texture", scene->textures,
            test_scn->textures, sel, update_list);
        draw_add_elem_glwidgets(scene, "mat", scene->materials,
            test_scn->materials, sel, update_list);
        draw_add_elem_glwidgets(
            scene, "shape", scene->shapes, test_scn->shapes, sel, update_list);
        draw_add_elem_glwidgets(scene, "instance", scene->instances,
            test_scn->instances, sel, update_list);
        draw_add_elem_glwidgets(
            scene, "node", scene->nodes, test_scn->nodes, sel, update_list);
        draw_add_elem_glwidgets(scene, "environment", scene->environments,
            test_scn->environments, sel, update_list);
        draw_add_elem_glwidgets(scene, "anim", scene->animations,
            test_scn->animations, sel, update_list);
    }
#endif
    return update_list.size() != update_len;
}

inline bool draw_glwidgets_scene_inspector(glwindow* win, const string& lbl,
    yocto_scene* scene, void*& sel, vector<tuple<string, void*>>& update_list,
    int height) {
    if (!scene || !sel) return false;
    begin_child_glwidget(win, "scrolling scene inspector", {0, height});

    auto update_len = update_list.size();

    for (auto camera : scene->cameras) {
        if (camera != sel) continue;
        if (draw_glwidgets_scene_inspector(win, camera, scene))
            update_list.push_back({"camera", camera});
    }
    for (auto shape : scene->shapes) {
        if (shape != sel) continue;
        if (draw_glwidgets_scene_inspector(win, shape, scene))
            update_list.push_back({"shape", shape});
    }
    for (auto texture : scene->textures) {
        if (texture != sel) continue;
        if (draw_glwidgets_scene_inspector(win, texture, scene))
            update_list.push_back({"texture", texture});
    }
    for (auto mat : scene->materials) {
        if (mat != sel) continue;
        if (draw_glwidgets_scene_inspector(win, mat, scene))
            update_list.push_back({"material", mat});
    }
    for (auto environment : scene->environments) {
        if (environment != sel) continue;
        if (draw_glwidgets_scene_inspector(win, environment, scene))
            update_list.push_back({"environment", environment});
    }
    for (auto instance : scene->instances) {
        if (instance != sel) continue;
        if (draw_glwidgets_scene_inspector(win, instance, scene))
            update_list.push_back({"instance", instance});
    }
    for (auto node : scene->nodes) {
        if (node != sel) continue;
        if (draw_glwidgets_scene_inspector(win, node, scene))
            update_list.push_back({"node", node});
    }
    for (auto animation : scene->animations) {
        if (animation != sel) continue;
        if (draw_glwidgets_scene_inspector(win, animation, scene))
            update_list.push_back({"animation", animation});
    }

    end_child_glwidget(win);
    return update_list.size() != update_len;
}

#endif
