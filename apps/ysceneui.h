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

inline const unordered_map<int, string>& animation_type_names() {
    static auto names = unordered_map<int, string>{
        {(int)yocto_interpolation_type::linear, "linear"},
        {(int)yocto_interpolation_type::step, "step"},
        {(int)yocto_interpolation_type::bezier, "bezier"},
    };
    return names;
}

template <typename T>
inline void draw_glwidgets_scene_tree(const glwindow& win, const string& lbl_,
    yocto_scene& scene, int index, const vector<T>& vals,
    tuple<string, int>& sel, const string& sel_type);

template <typename T>
inline void draw_glwidgets_scene_tree(const glwindow& win, const string& lbl_,
    yocto_scene& scene, int index, const vector<T*>& vals,
    tuple<string, int>& sel, const string& sel_type);

inline void draw_scene_tree_glwidgets_rec(const glwindow& win, const string& lbl_,
    yocto_scene& scene, const yocto_camera& value, tuple<string, int>& sel) {}
inline void draw_scene_tree_glwidgets_rec(const glwindow& win, const string& lbl_,
    yocto_scene& scene, const yocto_texture& value, tuple<string, int>& sel) {}
inline void draw_scene_tree_glwidgets_rec(const glwindow& win,
    const string& lbl_, yocto_scene& scene, const yocto_voltexture& value,
    tuple<string, int>& sel) {}
inline void draw_scene_tree_glwidgets_rec(const glwindow& win, const string& lbl_,
    yocto_scene& scene, const yocto_material& value, tuple<string, int>& sel) {
    draw_glwidgets_scene_tree(win, "emission", scene, value.emission_texture,
        scene.textures, sel, "texture");
    draw_glwidgets_scene_tree(win, "diffuse", scene, value.diffuse_texture,
        scene.textures, sel, "texture");
    draw_glwidgets_scene_tree(win, "specular", scene, value.specular_texture,
        scene.textures, sel, "texture");
    draw_glwidgets_scene_tree(
        win, "bump", scene, value.bump_texture, scene.textures, sel, "texture");
    draw_glwidgets_scene_tree(win, "displament", scene,
        value.displacement_texture, scene.textures, sel, "texture");
    draw_glwidgets_scene_tree(win, "normal", scene, value.normal_texture,
        scene.textures, sel, "texture");
}
inline void draw_scene_tree_glwidgets_rec(const glwindow& win, const string& lbl_,
    yocto_scene& scene, const yocto_shape& value, tuple<string, int>& sel) {
    draw_glwidgets_scene_tree(win, "material", scene, value.material,
        scene.materials, sel, "material");
}

inline void draw_scene_tree_glwidgets_rec(const glwindow& win, const string& lbl_,
    yocto_scene& scene, const yocto_instance& value, tuple<string, int>& sel) {
    draw_glwidgets_scene_tree(
        win, "shape", scene, value.shape, scene.shapes, sel, "shape");
}
inline void draw_scene_tree_glwidgets_rec(const glwindow& win,
    const string& lbl_, yocto_scene& scene, const yocto_environment& value,
    tuple<string, int>& sel) {
    draw_glwidgets_scene_tree(win, "emission", scene, value.emission_texture,
        scene.textures, sel, "texture");
}
inline void draw_scene_tree_glwidgets_rec(const glwindow& win, const string& lbl_,
    yocto_scene& scene, const yocto_scene_node& value, tuple<string, int>& sel) {
    draw_glwidgets_scene_tree(win, "instance", scene, value.instance,
        scene.instances, sel, "instance");
    draw_glwidgets_scene_tree(
        win, "camera", scene, value.camera, scene.cameras, sel, "camera");
    draw_glwidgets_scene_tree(win, "environment", scene, value.environment,
        scene.environments, sel, "environment");
    draw_glwidgets_scene_tree(
        win, "parent", scene, value.parent, scene.nodes, sel, "node");
    auto cid = 0;
    for (auto ch : value.children) {
        draw_glwidgets_scene_tree(win, "child" + to_string(cid++), scene, ch,
            scene.nodes, sel, "node");
    }
}
inline void draw_scene_tree_glwidgets_rec(const glwindow& win, const string& lbl_,
    yocto_scene& scene, const yocto_animation& value, tuple<string, int>& sel) {
    auto tid = 0;
    for (auto tg : value.node_targets) {
        draw_glwidgets_scene_tree(win, "target" + to_string(tid++), scene, tg,
            scene.nodes, sel, "node");
    }
}

template <typename T>
inline void draw_glwidgets_scene_tree(const glwindow& win, const string& lbl_,
    yocto_scene& scene, int index, const vector<T>& vals,
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
inline void draw_glwidgets_scene_tree(const glwindow& win, const string& lbl_,
    yocto_scene& scene, int index, const vector<T*>& vals,
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
    const glwindow& win, yocto_scene& scene, tuple<string, int>& sel) {
    if (!scene.cameras.empty() && begin_treenode_glwidget(win, "cameras")) {
        for (auto v = 0; v < scene.cameras.size(); v++)
            draw_glwidgets_scene_tree(
                win, "", scene, v, scene.cameras, sel, "camera");
        end_treenode_glwidget(win);
    }
    if (!scene.shapes.empty() && begin_treenode_glwidget(win, "shapes")) {
        for (auto v = 0; v < scene.shapes.size(); v++)
            draw_glwidgets_scene_tree(
                win, "", scene, v, scene.shapes, sel, "shape");
        end_treenode_glwidget(win);
    }
    if (!scene.instances.empty() && begin_treenode_glwidget(win, "instances")) {
        for (auto v = 0; v < scene.instances.size(); v++)
            draw_glwidgets_scene_tree(
                win, "", scene, v, scene.instances, sel, "instance");
        end_treenode_glwidget(win);
    }
    if (!scene.materials.empty() && begin_treenode_glwidget(win, "materials")) {
        for (auto v = 0; v < scene.materials.size(); v++)
            draw_glwidgets_scene_tree(
                win, "", scene, v, scene.materials, sel, "material");
        end_treenode_glwidget(win);
    }
    if (!scene.textures.empty() && begin_treenode_glwidget(win, "textures")) {
        for (auto v = 0; v < scene.textures.size(); v++)
            draw_glwidgets_scene_tree(
                win, "", scene, v, scene.textures, sel, "texture");
        end_treenode_glwidget(win);
    }
    if (!scene.environments.empty() &&
        begin_treenode_glwidget(win, "environments")) {
        for (auto v = 0; v < scene.environments.size(); v++)
            draw_glwidgets_scene_tree(
                win, "", scene, v, scene.environments, sel, "environment");
        end_treenode_glwidget(win);
    }
    if (!scene.nodes.empty() && begin_treenode_glwidget(win, "nodes")) {
        for (auto v = 0; v < scene.nodes.size(); v++)
            draw_glwidgets_scene_tree(
                win, "", scene, v, scene.nodes, sel, "node");
        end_treenode_glwidget(win);
    }
    if (!scene.animations.empty() && begin_treenode_glwidget(win, "animations")) {
        for (auto v = 0; v < scene.animations.size(); v++)
            draw_glwidgets_scene_tree(
                win, "", scene, v, scene.animations, sel, "animation");
        end_treenode_glwidget(win);
    }
}

/// Visit struct elements.
inline bool draw_glwidgets_scene_inspector(
    const glwindow& win, yocto_camera& value, yocto_scene& scene) {
    auto edited = 0;
    edited += draw_textinput_glwidget(win, "name", value.name);
    edited += draw_slider_glwidget(win, "frame.x", value.frame.x.x, -1, 1);
    edited += draw_slider_glwidget(win, "frame.y", value.frame.y.x, -1, 1);
    edited += draw_slider_glwidget(win, "frame.z", value.frame.z.x, -1, 1);
    edited += draw_slider_glwidget(win, "frame.o", value.frame.o.x, -10, 10);
    edited += draw_checkbox_glwidget(win, "ortho", value.orthographic);
    edited += draw_slider_glwidget(win, "film", value.film_size, 0.01f, 1);
    edited += draw_slider_glwidget(win, "focal", value.focal_length, 0.01f, 1);
    edited += draw_slider_glwidget(
        win, "focus", value.focus_distance, 0.01f, 1000);
    edited += draw_slider_glwidget(win, "aperture", value.lens_aperture, 0, 5);
    return edited;
}

/// Visit struct elements.
inline bool draw_glwidgets_scene_inspector(
    const glwindow& win, yocto_texture& value, yocto_scene& scene) {
    auto edited = 0;
    edited += draw_textinput_glwidget(win, "name", value.name);
    edited += draw_textinput_glwidget(win, "path", value.filename);
    edited += draw_checkbox_glwidget(win, "clamp_to_edge", value.clamp_to_edge);
    edited += draw_slider_glwidget(win, "scale", value.height_scale, 0, 1);
    edited += draw_checkbox_glwidget(win, "ldr_as_linear", value.ldr_as_linear);
    draw_label_glwidgets(win, "hdr_image", "%d x %d", value.hdr_image.width,
        value.hdr_image.height);
    draw_label_glwidgets(win, "ldr_image", "%d x %d", value.ldr_image.width,
        value.ldr_image.height);
    return edited;
}

inline bool draw_glwidgets_scene_inspector(
    const glwindow& win, yocto_material& value, yocto_scene& scene) {
    auto edited = 0;
    edited += draw_textinput_glwidget(win, "name", value.name);
    edited += draw_coloredit_glwidget(win, "emission", value.emission);  // TODO:
                                                                         // HDR
    edited += draw_coloredit_glwidget(win, "diffuse", value.diffuse);
    edited += draw_coloredit_glwidget(win, "specular", value.specular);
    edited += draw_coloredit_glwidget(win, "transmission", value.transmission);
    edited += draw_slider_glwidget(win, "roughness", value.roughness, 0, 1);
    edited += draw_slider_glwidget(win, "opacity", value.opacity, 0, 1);
    edited += draw_checkbox_glwidget(win, "fresnel", value.fresnel);
    continue_glwidgets_line(win);
    edited += draw_checkbox_glwidget(win, "refract", value.refract);
    edited += draw_coloredit_glwidget(
        win, "volume_density", value.volume_density);  // 0, 10
    edited += draw_coloredit_glwidget(
        win, "volume_albedo", value.volume_albedo);  // 0, 1
    edited += draw_slider_glwidget(
        win, "volume_phaseg", value.volume_phaseg, -1, 1);
    edited += draw_combobox_glwidget(
        win, "emission_texture", value.emission_texture, scene.textures, true);
    edited += draw_combobox_glwidget(
        win, "diffuse_texture", value.diffuse_texture, scene.textures, true);
    edited += draw_combobox_glwidget(
        win, "specular_texture", value.specular_texture, scene.textures, true);
    edited += draw_combobox_glwidget(win, "transmission_texture",
        value.transmission_texture, scene.textures, true);
    edited += draw_combobox_glwidget(
        win, "opacity_texture", value.opacity_texture, scene.textures, true);
    edited += draw_combobox_glwidget(win, "roughness_texture",
        value.roughness_texture, scene.textures, true);
    edited += draw_combobox_glwidget(
        win, "bump_texture", value.bump_texture, scene.textures, true);
    edited += draw_combobox_glwidget(win, "displacement_texture",
        value.displacement_texture, scene.textures, true);
    edited += draw_combobox_glwidget(
        win, "normal_texture", value.normal_texture, scene.textures, true);
    edited += draw_checkbox_glwidget(win, "base metallic", value.base_metallic);
    edited += draw_checkbox_glwidget(win, "glTF textures", value.gltf_textures);
    return edited;
}

inline bool draw_glwidgets_scene_inspector(
    const glwindow& win, yocto_shape& value, yocto_scene& scene) {
    auto edited = 0;
    edited += draw_textinput_glwidget(win, "name", value.name);
    edited += draw_textinput_glwidget(win, "path", value.filename);
    edited += draw_combobox_glwidget(
        win, "material", value.material, scene.materials, true);
    draw_label_glwidgets(win, "lines", "%ld", value.lines.size());
    draw_label_glwidgets(win, "triangles", "%ld", value.triangles.size());
    draw_label_glwidgets(win, "quads", "%ld", value.quads.size());
    draw_label_glwidgets(win, "pos", "%ld", value.positions.size());
    draw_label_glwidgets(win, "norm", "%ld", value.normals.size());
    draw_label_glwidgets(win, "texcoord", "%ld", value.texturecoords.size());
    draw_label_glwidgets(win, "color", "%ld", value.colors.size());
    draw_label_glwidgets(win, "radius", "%ld", value.radius.size());
    draw_label_glwidgets(win, "tangsp", "%ld", value.tangentspaces.size());
    return edited;
}

inline bool draw_glwidgets_scene_inspector(
    const glwindow& win, yocto_instance& value, yocto_scene& scene) {
    auto edited = 0;
    edited += draw_textinput_glwidget(win, "name", value.name);
    edited += draw_slider_glwidget(win, "frame.x", value.frame.x, -1, 1);
    edited += draw_slider_glwidget(win, "frame.y", value.frame.y, -1, 1);
    edited += draw_slider_glwidget(win, "frame.z", value.frame.z, -1, 1);
    edited += draw_slider_glwidget(win, "frame.o", value.frame.o, -10, 10);
    edited += draw_combobox_glwidget(
        win, "shape", value.shape, scene.shapes, true);
    return edited;
}

inline bool draw_glwidgets_scene_inspector(
    const glwindow& win, yocto_environment& value, yocto_scene& scene) {
    auto edited = 0;
    edited += draw_textinput_glwidget(win, "name", value.name);
    edited += draw_slider_glwidget(win, "frame.x", value.frame.x, -1, 1);
    edited += draw_slider_glwidget(win, "frame.y", value.frame.y, -1, 1);
    edited += draw_slider_glwidget(win, "frame.z", value.frame.z, -1, 1);
    edited += draw_slider_glwidget(win, "frame.o", value.frame.o, -10, 10);
    edited += draw_coloredit_glwidget(win, "ke", value.emission);  // TODO: HDR
    edited += draw_combobox_glwidget(
        win, "ke texture", value.emission_texture, scene.textures, true);
    return edited;
}

inline bool draw_glwidgets_scene_inspector(
    const glwindow& win, yocto_scene_node& value, yocto_scene& scene) {
    auto edited = 0;
    edited += draw_textinput_glwidget(win, "name", value.name);
    edited += draw_combobox_glwidget(
        win, "parent", value.parent, scene.nodes, true);
    edited += draw_slider_glwidget(win, "local.x", value.local.x, -1, 1);
    edited += draw_slider_glwidget(win, "local.y", value.local.y, -1, 1);
    edited += draw_slider_glwidget(win, "local.z", value.local.z, -1, 1);
    edited += draw_slider_glwidget(win, "local.o", value.local.o, -10, 10);
    edited += draw_slider_glwidget(
        win, "translation", value.translation, -10, 10);
    edited += draw_slider_glwidget(win, "rotation", value.rotation, -1, 1);
    edited += draw_slider_glwidget(win, "scale", value.scale, 0, 10);
    edited += draw_combobox_glwidget(
        win, "camera", value.camera, scene.cameras, true);
    edited += draw_combobox_glwidget(
        win, "instance", value.instance, scene.instances, true);
    edited += draw_combobox_glwidget(
        win, "environment", value.environment, scene.environments, true);
    return edited;
}

inline bool draw_glwidgets_scene_inspector(
    const glwindow& win, yocto_animation& value, yocto_scene& scene) {
    auto edited = 0;
    edited += draw_textinput_glwidget(win, "name", value.name);
    edited += draw_textinput_glwidget(win, "path", value.filename);
    edited += draw_textinput_glwidget(win, "group", value.animation_group);
    // edited += draw_combobox_glwidget(win, "type", &value.type,
    // animation_type_names());
    draw_label_glwidgets(win, "times", "%ld", value.keyframes_times.size());
    draw_label_glwidgets(
        win, "translation", "%ld", value.translation_keyframes.size());
    draw_label_glwidgets(win, "rotation", "%ld", value.rotation_keyframes.size());
    draw_label_glwidgets(win, "scale", "%ld", value.scale_keyframes.size());
    draw_label_glwidgets(
        win, "weights", "%ld", value.morph_weights_keyframes.size());
    draw_label_glwidgets(win, "targets", "%ld", value.node_targets.size());
    return edited;
}

inline bool draw_glwidgets_scene_tree(const glwindow& win, const string& lbl,
    yocto_scene& scene, tuple<string, int>& sel,
    vector<tuple<string, int>>& update_list, int height) {
    draw_glwidgets_scene_tree(win, scene, sel);
    auto update_len = update_list.size();
#if 0
    if (test_scn) {
        draw_add_elem_glwidgets(
            scene, "camera", scene.cameras, test_scn->cameras, sel, update_list);
        draw_add_elem_glwidgets(scene, "texture", scene.textures,
            test_scn->textures, sel, update_list);
        draw_add_elem_glwidgets(scene, "mat", scene.materials,
            test_scn->materials, sel, update_list);
        draw_add_elem_glwidgets(
            scene, "shape", scene.shapes, test_scn->shapes, sel, update_list);
        draw_add_elem_glwidgets(scene, "instance", scene.instances,
            test_scn->instances, sel, update_list);
        draw_add_elem_glwidgets(
            scene, "node", scene.nodes, test_scn->nodes, sel, update_list);
        draw_add_elem_glwidgets(scene, "environment", scene.environments,
            test_scn->environments, sel, update_list);
        draw_add_elem_glwidgets(scene, "anim", scene.animations,
            test_scn->animations, sel, update_list);
    }
#endif
    return update_list.size() != update_len;
}

inline bool draw_glwidgets_scene_inspector(const glwindow& win,
    const string& lbl, yocto_scene& scene, tuple<string, int>& sel,
    vector<tuple<string, int>>& update_list, int height) {
    if (get<0>(sel) == "") return false;
    begin_child_glwidget(win, "scrolling scene inspector", {0, height});

    auto update_len = update_list.size();

    if (get<0>(sel) == "camera")
        if (draw_glwidgets_scene_inspector(win, scene.cameras[get<1>(sel)], scene))
            update_list.push_back({"camera", get<1>(sel)});
    if (get<0>(sel) == "shape")
        if (draw_glwidgets_scene_inspector(win, scene.shapes[get<1>(sel)], scene))
            update_list.push_back({"shape", get<1>(sel)});
    if (get<0>(sel) == "texture")
        if (draw_glwidgets_scene_inspector(
                win, scene.textures[get<1>(sel)], scene))
            update_list.push_back({"texture", get<1>(sel)});
    if (get<0>(sel) == "material")
        if (draw_glwidgets_scene_inspector(
                win, scene.materials[get<1>(sel)], scene))
            update_list.push_back({"material", get<1>(sel)});
    if (get<0>(sel) == "environment")
        if (draw_glwidgets_scene_inspector(
                win, scene.environments[get<1>(sel)], scene))
            update_list.push_back({"environment", get<1>(sel)});
    if (get<0>(sel) == "instance")
        if (draw_glwidgets_scene_inspector(
                win, scene.instances[get<1>(sel)], scene))
            update_list.push_back({"instance", get<1>(sel)});
    if (get<0>(sel) == "node")
        if (draw_glwidgets_scene_inspector(win, scene.nodes[get<1>(sel)], scene))
            update_list.push_back({"node", get<1>(sel)});
    if (get<0>(sel) == "animation")
        if (draw_glwidgets_scene_inspector(
                win, scene.animations[get<1>(sel)], scene))
            update_list.push_back({"animation", get<1>(sel)});

    end_child_glwidget(win);
    return update_list.size() != update_len;
}

#endif
