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
inline void draw_glwidgets_scene_tree(glwindow* win, const string& lbl_,
    yocto_scene* scene, int index, const vector<T>& vals,
                                      tuple<string, int>& sel, const string& sel_type);

template <typename T>
inline void draw_glwidgets_scene_tree(glwindow* win, const string& lbl_,
    yocto_scene* scene, int index, const vector<T*>& vals,
                                      tuple<string, int>& sel, const string& sel_type);

inline void draw_scene_tree_glwidgets_rec(glwindow* win, const string& lbl_,
    yocto_scene* scene, const yocto_camera& val, tuple<string, int>& sel) {}
inline void draw_scene_tree_glwidgets_rec(glwindow* win, const string& lbl_,
    yocto_scene* scene, yocto_texture* val, tuple<string, int>& sel) {}
inline void draw_scene_tree_glwidgets_rec(glwindow* win, const string& lbl_,
    yocto_scene* scene, yocto_voltexture* val, tuple<string, int>& sel) {}
inline void draw_scene_tree_glwidgets_rec(glwindow* win, const string& lbl_,
    yocto_scene* scene, yocto_material* val, tuple<string, int>& sel) {
    draw_glwidgets_scene_tree(win, "emission", scene, val->emission_texture,
        scene->textures, sel, "texture");
    draw_glwidgets_scene_tree(win, "diffuse", scene, val->diffuse_texture,
        scene->textures, sel, "texture");
    draw_glwidgets_scene_tree(win, "specular", scene, val->specular_texture,
        scene->textures, sel, "texture");
    draw_glwidgets_scene_tree(
        win, "bump", scene, val->bump_texture, scene->textures, sel, "texture");
    draw_glwidgets_scene_tree(win, "displament", scene,
        val->displacement_texture, scene->textures, sel, "texture");
    draw_glwidgets_scene_tree(win, "normal", scene, val->normal_texture,
        scene->textures, sel, "texture");
}
inline void draw_scene_tree_glwidgets_rec(glwindow* win, const string& lbl_,
    yocto_scene* scene, yocto_shape* val, tuple<string, int>& sel) {
    draw_glwidgets_scene_tree(win, "material", scene, val->material,
        scene->materials, sel, "material");
}

inline void draw_scene_tree_glwidgets_rec(glwindow* win, const string& lbl_,
    yocto_scene* scene, yocto_instance* val, tuple<string, int>& sel) {
    draw_glwidgets_scene_tree(
        win, "shape", scene, val->shape, scene->shapes, sel, "shape");
}
inline void draw_scene_tree_glwidgets_rec(glwindow* win, const string& lbl_,
    yocto_scene* scene, yocto_environment* val, tuple<string, int>& sel) {
    draw_glwidgets_scene_tree(win, "emission", scene, val->emission_texture,
        scene->textures, sel, "texture");
}
inline void draw_scene_tree_glwidgets_rec(glwindow* win, const string& lbl_,
    yocto_scene* scene, yocto_scene_node* val, tuple<string, int>& sel) {
    draw_glwidgets_scene_tree(win, "instance", scene, val->instance,
        scene->instances, sel, "instance");
    draw_glwidgets_scene_tree(
        win, "camera", scene, val->camera, scene->cameras_, sel, "camera");
    draw_glwidgets_scene_tree(win, "environment", scene, val->environment,
        scene->environments, sel, "environment");
    draw_glwidgets_scene_tree(
        win, "parent", scene, val->parent, scene->nodes, sel, "node");
    auto cid = 0;
    for (auto ch : val->children) {
        draw_glwidgets_scene_tree(win, "child" + to_string(cid++), scene, ch,
            scene->nodes, sel, "node");
    }
}
inline void draw_scene_tree_glwidgets_rec(glwindow* win, const string& lbl_,
    yocto_scene* scene, yocto_animation* val, tuple<string, int>& sel) {
    auto tid = 0;
    for (auto tg : val->node_targets) {
        draw_glwidgets_scene_tree(win, "target" + to_string(tid++), scene, tg,
            scene->nodes, sel, "node");
    }
}

template <typename T>
inline void draw_glwidgets_scene_tree(glwindow* win, const string& lbl_,
                                      yocto_scene* scene, int index, const vector<T>& vals,
                                      tuple<string, int>& sel, const string& sel_type) {
    if (index < 0) return;
    auto lbl = vals[index].name;
    if (!lbl_.empty()) lbl = lbl_ + ": " + vals[index].name;
    auto selected = sel == tuple<string, int>{sel_type, index};
    if (begin_selectabletreenode_glwidget(win, lbl.c_str(), selected)) {
        draw_scene_tree_glwidgets_rec(win, lbl_, scene, vals[index], sel);
        end_treenode_glwidget(win);
    }
    if (selected) sel = {sel_type, index};
}

template <typename T>
inline void draw_glwidgets_scene_tree(glwindow* win, const string& lbl_,
                                      yocto_scene* scene, int index, const vector<T*>& vals,
                                      tuple<string, int>& sel, const string& sel_type) {
    if (index < 0) return;
    auto lbl = vals[index]->name;
    if (!lbl_.empty()) lbl = lbl_ + ": " + vals[index]->name;
    auto selected = sel == tuple<string, int>{sel_type, index};
    if (begin_selectabletreenode_glwidget(win, lbl.c_str(), selected)) {
        draw_scene_tree_glwidgets_rec(win, lbl_, scene, vals[index], sel);
        end_treenode_glwidget(win);
    }
    if (selected) sel = {sel_type, index};
}

inline void draw_glwidgets_scene_tree(
    glwindow* win, yocto_scene* scene, tuple<string, int>& sel) {
    if (!scene->cameras_.empty() && begin_treenode_glwidget(win, "cameras")) {
        for (auto v = 0; v < scene->cameras_.size(); v++)
            draw_glwidgets_scene_tree(
                win, "", scene, v, scene->cameras_, sel, "camera");
        end_treenode_glwidget(win);
    }
    if (!scene->shapes.empty() && begin_treenode_glwidget(win, "shapes")) {
        for (auto v = 0; v < scene->shapes.size(); v++)
            draw_glwidgets_scene_tree(
                win, "", scene, v, scene->shapes, sel, "shape");
        end_treenode_glwidget(win);
    }
    if (!scene->instances.empty() && begin_treenode_glwidget(win, "instances")) {
        for (auto v = 0; v < scene->instances.size(); v++)
            draw_glwidgets_scene_tree(
                win, "", scene, v, scene->instances, sel, "instance");
        end_treenode_glwidget(win);
    }
    if (!scene->materials.empty() && begin_treenode_glwidget(win, "materials")) {
        for (auto v = 0; v < scene->materials.size(); v++)
            draw_glwidgets_scene_tree(
                win, "", scene, v, scene->materials, sel, "material");
        end_treenode_glwidget(win);
    }
    if (!scene->textures.empty() && begin_treenode_glwidget(win, "textures")) {
        for (auto v = 0; v < scene->textures.size(); v++)
            draw_glwidgets_scene_tree(
                win, "", scene, v, scene->textures, sel, "texture");
        end_treenode_glwidget(win);
    }
    if (!scene->environments.empty() &&
        begin_treenode_glwidget(win, "environments")) {
        for (auto v = 0; v < scene->environments.size(); v++)
            draw_glwidgets_scene_tree(
                win, "", scene, v, scene->environments, sel, "environment");
        end_treenode_glwidget(win);
    }
    if (!scene->nodes.empty() && begin_treenode_glwidget(win, "nodes")) {
        for (auto v = 0; v < scene->nodes.size(); v++)
            draw_glwidgets_scene_tree(
                win, "", scene, v, scene->nodes, sel, "node");
        end_treenode_glwidget(win);
    }
    if (!scene->animations.empty() && begin_treenode_glwidget(win, "animations")) {
        for (auto v = 0; v < scene->animations.size(); v++)
            draw_glwidgets_scene_tree(
                win, "", scene, v, scene->animations, sel, "animation");
        end_treenode_glwidget(win);
    }
}

/// Visit struct elements.
inline bool draw_glwidgets_scene_inspector(
    glwindow* win, yocto_camera& val, yocto_scene* scene) {
    auto edited = 0;
    edited += draw_textinput_glwidget(win, "name", val.name);
    edited += draw_slider_glwidget(win, "frame.x", val.frame.x.x, -1, 1);
    edited += draw_slider_glwidget(win, "frame.y", val.frame.y.x, -1, 1);
    edited += draw_slider_glwidget(win, "frame.z", val.frame.z.x, -1, 1);
    edited += draw_slider_glwidget(win, "frame.o", val.frame.o.x, -10, 10);
    edited += draw_checkbox_glwidget(win, "ortho", val.orthographic);
    edited += draw_slider_glwidget(win, "film", val.film_size, 0.01f, 1);
    edited += draw_slider_glwidget(win, "focal", val.focal_length, 0.01f, 1);
    edited += draw_slider_glwidget(win, "focus", val.focus_distance, 0.01f, 1000);
    edited += draw_slider_glwidget(win, "aperture", val.lens_aperture, 0, 5);
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
        win, "camera", val->camera, scene->cameras_, true);
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
    yocto_scene* scene, tuple<string, int>& sel,
    vector<tuple<string, int>>& update_list, int height) {
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
    yocto_scene* scene, tuple<string, int>& sel,
    vector<tuple<string, int>>& update_list, int height) {
    if (!scene || get<0>(sel) == "") return false;
    begin_child_glwidget(win, "scrolling scene inspector", {0, height});

    auto update_len = update_list.size();

    if (get<0>(sel) == "camera")
        if (draw_glwidgets_scene_inspector(
                win, scene->cameras_[get<1>(sel)], scene))
            update_list.push_back({"camera", get<1>(sel)});
    if (get<0>(sel) == "shape")
        if (draw_glwidgets_scene_inspector(win, scene->shapes[get<1>(sel)], scene))
            update_list.push_back({"shape", get<1>(sel)});
    if (get<0>(sel) == "texture")
        if (draw_glwidgets_scene_inspector(
                win, scene->textures[get<1>(sel)], scene))
            update_list.push_back({"texture", get<1>(sel)});
    if (get<0>(sel) == "material")
        if (draw_glwidgets_scene_inspector(
                win, scene->materials[get<1>(sel)], scene))
            update_list.push_back({"material", get<1>(sel)});
    if (get<0>(sel) == "environment")
        if (draw_glwidgets_scene_inspector(
                win, scene->environments[get<1>(sel)], scene))
            update_list.push_back({"environment", get<1>(sel)});
    if (get<0>(sel) == "instance")
        if (draw_glwidgets_scene_inspector(
                win, scene->instances[get<1>(sel)], scene))
            update_list.push_back({"instance", get<1>(sel)});
    if (get<0>(sel) == "node")
        if (draw_glwidgets_scene_inspector(win, scene->nodes[get<1>(sel)], scene))
            update_list.push_back({"node", get<1>(sel)});
    if (get<0>(sel) == "animation")
        if (draw_glwidgets_scene_inspector(
                win, scene->animations[get<1>(sel)], scene))
            update_list.push_back({"animation", get<1>(sel)});

    end_child_glwidget(win);
    return update_list.size() != update_len;
}

#endif
