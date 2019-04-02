//
// Utilities to display a scene graph using ImGui.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2019 Fabio Pellacini
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

#ifndef _YOCTO_SCENEUI_H_
#define _YOCTO_SCENEUI_H_

#include "../yocto/yocto_scene.h"
#include "yocto_opengl.h"
using namespace yocto;

#include <typeindex>
#include <unordered_map>
using std::type_index;
using std::unordered_map;

inline const unordered_map<int, string>& animation_type_names() {
    static auto names = unordered_map<int, string>{
        {(int)yocto_interpolation_type::linear, "linear"},
        {(int)yocto_interpolation_type::step, "step"},
        {(int)yocto_interpolation_type::bezier, "bezier"},
    };
    return names;
}

template <typename T>
inline void draw_opengl_widgets_scene_tree(const opengl_window& win,
    const string& lbl_, yocto_scene& scene, int index, const vector<T>& vals,
    pair<type_index, int>& sel);

template <typename T>
inline void draw_opengl_widgets_scene_tree(const opengl_window& win,
    const string& lbl_, yocto_scene& scene, int index, const vector<T*>& vals,
    pair<type_index, int>& sel);

inline void draw_scene_tree_opengl_widgets_rec(const opengl_window& win,
    const string& lbl_, yocto_scene& scene, const yocto_camera& value,
    pair<type_index, int>& sel) {}
inline void draw_scene_tree_opengl_widgets_rec(const opengl_window& win,
    const string& lbl_, yocto_scene& scene, const yocto_texture& value,
    pair<type_index, int>& sel) {}
inline void draw_scene_tree_opengl_widgets_rec(const opengl_window& win,
    const string& lbl_, yocto_scene& scene, const yocto_voltexture& value,
    pair<type_index, int>& sel) {}
inline void draw_scene_tree_opengl_widgets_rec(const opengl_window& win,
    const string& lbl_, yocto_scene& scene, const yocto_material& value,
    pair<type_index, int>& sel) {
    draw_opengl_widgets_scene_tree(
        win, "emission", scene, value.emission_texture, scene.textures, sel);
    draw_opengl_widgets_scene_tree(
        win, "diffuse", scene, value.diffuse_texture, scene.textures, sel);
    draw_opengl_widgets_scene_tree(
        win, "specular", scene, value.specular_texture, scene.textures, sel);
    draw_opengl_widgets_scene_tree(win, "displament", scene,
        value.displacement_texture, scene.textures, sel);
    draw_opengl_widgets_scene_tree(
        win, "normal", scene, value.normal_texture, scene.textures, sel);
}
inline void draw_scene_tree_opengl_widgets_rec(const opengl_window& win,
    const string& lbl_, yocto_scene& scene, const yocto_shape& value,
    pair<type_index, int>& sel) {}

inline void draw_scene_tree_opengl_widgets_rec(const opengl_window& win,
    const string& lbl_, yocto_scene& scene, const yocto_instance& value,
    pair<type_index, int>& sel) {
    draw_opengl_widgets_scene_tree(
        win, "shape", scene, value.shape, scene.shapes, sel);
    draw_opengl_widgets_scene_tree(
        win, "material", scene, value.material, scene.materials, sel);
}
inline void draw_scene_tree_opengl_widgets_rec(const opengl_window& win,
    const string& lbl_, yocto_scene& scene, const yocto_environment& value,
    pair<type_index, int>& sel) {
    draw_opengl_widgets_scene_tree(
        win, "emission", scene, value.emission_texture, scene.textures, sel);
}
inline void draw_scene_tree_opengl_widgets_rec(const opengl_window& win,
    const string& lbl_, yocto_scene& scene, const yocto_scene_node& value,
    pair<type_index, int>& sel) {
    draw_opengl_widgets_scene_tree(
        win, "instance", scene, value.instance, scene.instances, sel);
    draw_opengl_widgets_scene_tree(
        win, "camera", scene, value.camera, scene.cameras, sel);
    draw_opengl_widgets_scene_tree(
        win, "environment", scene, value.environment, scene.environments, sel);
    draw_opengl_widgets_scene_tree(
        win, "parent", scene, value.parent, scene.nodes, sel);
    auto cid = 0;
    for (auto ch : value.children) {
        draw_opengl_widgets_scene_tree(
            win, "child" + std::to_string(cid++), scene, ch, scene.nodes, sel);
    }
}
inline void draw_scene_tree_opengl_widgets_rec(const opengl_window& win,
    const string& lbl_, yocto_scene& scene, const yocto_animation& value,
    pair<type_index, int>& sel) {
    auto tid = 0;
    for (auto tg : value.node_targets) {
        draw_opengl_widgets_scene_tree(
            win, "target" + std::to_string(tid++), scene, tg, scene.nodes, sel);
    }
}

template <typename T>
inline void draw_opengl_widgets_scene_tree(const opengl_window& win,
    const string& lbl_, yocto_scene& scene, int index, const vector<T>& vals,
    pair<type_index, int>& sel) {
    if (index < 0) return;
    auto lbl = vals[index].name;
    if (!empty(lbl_)) lbl = lbl_ + ": " + vals[index].name;
    auto selected = sel == pair<type_index, int>{type_index(typeid(T)), index};
    if (begin_selectabletreenode_opengl_widget(win, lbl.c_str(), selected)) {
        draw_scene_tree_opengl_widgets_rec(win, lbl_, scene, vals[index], sel);
        end_treenode_opengl_widget(win);
    }
    if (selected) sel = {type_index(typeid(T)), index};
}

template <typename T>
inline void draw_opengl_widgets_scene_tree(const opengl_window& win,
    const string& lbl_, yocto_scene& scene, int index, const vector<T*>& vals,
    pair<type_index, int>& sel) {
    if (index < 0) return;
    auto lbl = vals[index]->name;
    if (!empty(lbl_)) lbl = lbl_ + ": " + vals[index]->name;
    auto selected = sel == pair<type_index, int>{type_index(typeid(T)), index};
    if (begin_selectabletreenode_opengl_widget(win, lbl.c_str(), selected)) {
        draw_scene_tree_opengl_widgets_rec(win, lbl_, scene, vals[index], sel);
        end_treenode_opengl_widget(win);
    }
    if (selected) sel = {type_index(typeid(T)), index};
}

inline void draw_opengl_widgets_scene_tree(
    const opengl_window& win, yocto_scene& scene, pair<type_index, int>& sel) {
    if (!empty(scene.cameras) && begin_treenode_opengl_widget(win, "cameras")) {
        for (auto v = 0; v < scene.cameras.size(); v++)
            draw_opengl_widgets_scene_tree(
                win, "", scene, v, scene.cameras, sel);
        end_treenode_opengl_widget(win);
    }
    if (!empty(scene.shapes) && begin_treenode_opengl_widget(win, "shapes")) {
        for (auto v = 0; v < scene.shapes.size(); v++)
            draw_opengl_widgets_scene_tree(
                win, "", scene, v, scene.shapes, sel);
        end_treenode_opengl_widget(win);
    }
    if (!empty(scene.instances) &&
        begin_treenode_opengl_widget(win, "instances")) {
        for (auto v = 0; v < scene.instances.size(); v++)
            draw_opengl_widgets_scene_tree(
                win, "", scene, v, scene.instances, sel);
        end_treenode_opengl_widget(win);
    }
    if (!empty(scene.materials) &&
        begin_treenode_opengl_widget(win, "materials")) {
        for (auto v = 0; v < scene.materials.size(); v++)
            draw_opengl_widgets_scene_tree(
                win, "", scene, v, scene.materials, sel);
        end_treenode_opengl_widget(win);
    }
    if (!empty(scene.textures) &&
        begin_treenode_opengl_widget(win, "textures")) {
        for (auto v = 0; v < scene.textures.size(); v++)
            draw_opengl_widgets_scene_tree(
                win, "", scene, v, scene.textures, sel);
        end_treenode_opengl_widget(win);
    }
    if (!empty(scene.environments) &&
        begin_treenode_opengl_widget(win, "environments")) {
        for (auto v = 0; v < scene.environments.size(); v++)
            draw_opengl_widgets_scene_tree(
                win, "", scene, v, scene.environments, sel);
        end_treenode_opengl_widget(win);
    }
    if (!empty(scene.nodes) && begin_treenode_opengl_widget(win, "nodes")) {
        for (auto v = 0; v < scene.nodes.size(); v++)
            draw_opengl_widgets_scene_tree(win, "", scene, v, scene.nodes, sel);
        end_treenode_opengl_widget(win);
    }
    if (!empty(scene.animations) &&
        begin_treenode_opengl_widget(win, "animations")) {
        for (auto v = 0; v < scene.animations.size(); v++)
            draw_opengl_widgets_scene_tree(
                win, "", scene, v, scene.animations, sel);
        end_treenode_opengl_widget(win);
    }
}

/// Visit struct elements.
inline bool draw_opengl_widgets_scene_inspector(
    const opengl_window& win, yocto_camera& value, yocto_scene& scene) {
    auto edited = 0;
    edited += draw_textinput_opengl_widget(win, "name", value.name);
    edited += draw_slider_opengl_widget(win, "frame.x", value.frame.x, -1, 1);
    edited += draw_slider_opengl_widget(win, "frame.y", value.frame.y, -1, 1);
    edited += draw_slider_opengl_widget(win, "frame.z", value.frame.z, -1, 1);
    edited += draw_slider_opengl_widget(win, "frame.o", value.frame.o, -10, 10);
    edited += draw_checkbox_opengl_widget(win, "ortho", value.orthographic);
    edited += draw_slider_opengl_widget(
        win, "film width", value.film_width, 0.01f, 1);
    edited += draw_slider_opengl_widget(
        win, "film height", value.film_height, 0.01f, 1);
    edited += draw_slider_opengl_widget(
        win, "focal length", value.focal_length, 0.01f, 1);
    edited += draw_slider_opengl_widget(
        win, "focus distance", value.focus_distance, 0.01f, 1000);
    edited += draw_slider_opengl_widget(
        win, "lens aperture", value.lens_aperture, 0, 5);
    auto from = value.frame.o,
         to   = value.frame.o - value.focus_distance * value.frame.z;
    draw_slider_opengl_widget(win, "!!from", from, -10, 10);
    draw_slider_opengl_widget(win, "!!to", to, -10, 10);
    return edited;
}

/// Visit struct elements.
inline bool draw_opengl_widgets_scene_inspector(
    const opengl_window& win, yocto_texture& value, yocto_scene& scene) {
    auto edited = 0;
    edited += draw_textinput_opengl_widget(win, "name", value.name);
    edited += draw_textinput_opengl_widget(win, "path", value.filename);
    edited += draw_checkbox_opengl_widget(
        win, "clamp_to_edge", value.clamp_to_edge);
    edited += draw_slider_opengl_widget(win, "scale", value.height_scale, 0, 1);
    edited += draw_checkbox_opengl_widget(
        win, "ldr_as_linear", value.ldr_as_linear);
    draw_label_opengl_widget(win, "hdr_image", "%d x %d",
        value.hdr_image.size().x, value.hdr_image.size().y);
    draw_label_opengl_widget(win, "ldr_image", "%d x %d",
        value.ldr_image.size().x, value.ldr_image.size().y);
    return edited;
}

inline bool draw_opengl_widgets_scene_inspector(
    const opengl_window& win, yocto_material& value, yocto_scene& scene) {
    auto edited = 0;
    edited += draw_textinput_opengl_widget(win, "name", value.name);
    edited += draw_hdr_coloredit_opengl_widget(win, "emission", value.emission);
    edited += draw_coloredit_opengl_widget(win, "diffuse", value.diffuse);
    edited += draw_coloredit_opengl_widget(win, "specular", value.specular);
    edited += draw_coloredit_opengl_widget(
        win, "transmission", value.transmission);
    edited += draw_slider_opengl_widget(
        win, "roughness", value.roughness, 0, 1);
    edited += draw_slider_opengl_widget(win, "opacity", value.opacity, 0, 1);
    continue_opengl_widget_line(win);
    edited += draw_checkbox_opengl_widget(win, "fresnel", value.fresnel);
    continue_opengl_widget_line(win);
    edited += draw_checkbox_opengl_widget(win, "refract", value.refract);
    edited += draw_coloredit_opengl_widget(win, "volume_density",
        value.volume_density);  // 0, 10
    edited += draw_coloredit_opengl_widget(win, "volume_albedo",
        value.volume_albedo);  // 0, 1
    edited += draw_slider_opengl_widget(
        win, "volume_phaseg", value.volume_phaseg, -1, 1);
    edited += draw_combobox_opengl_widget(
        win, "emission_texture", value.emission_texture, scene.textures, true);
    edited += draw_combobox_opengl_widget(
        win, "diffuse_texture", value.diffuse_texture, scene.textures, true);
    edited += draw_combobox_opengl_widget(
        win, "specular_texture", value.specular_texture, scene.textures, true);
    edited += draw_combobox_opengl_widget(win, "transmission_texture",
        value.transmission_texture, scene.textures, true);
    edited += draw_combobox_opengl_widget(win, "roughness_texture",
        value.roughness_texture, scene.textures, true);
    edited += draw_combobox_opengl_widget(win, "displacement_texture",
        value.displacement_texture, scene.textures, true);
    edited += draw_combobox_opengl_widget(
        win, "normal_texture", value.normal_texture, scene.textures, true);
    edited += draw_checkbox_opengl_widget(
        win, "base metallic", value.base_metallic);
    edited += draw_checkbox_opengl_widget(
        win, "glTF textures", value.gltf_textures);
    return edited;
}

inline bool draw_opengl_widgets_scene_inspector(
    const opengl_window& win, yocto_shape& value, yocto_scene& scene) {
    auto edited = 0;
    edited += draw_textinput_opengl_widget(win, "name", value.name);
    edited += draw_textinput_opengl_widget(win, "path", value.filename);
    draw_label_opengl_widget(win, "lines", "%ld", value.lines.size());
    draw_label_opengl_widget(win, "triangles", "%ld", value.triangles.size());
    draw_label_opengl_widget(win, "quads", "%ld", value.quads.size());
    draw_label_opengl_widget(
        win, "quads pos", "%ld", value.quads_positions.size());
    draw_label_opengl_widget(
        win, "quads norm", "%ld", value.quads_normals.size());
    draw_label_opengl_widget(
        win, "quads texcoord", "%ld", value.quads_texturecoords.size());
    draw_label_opengl_widget(win, "pos", "%ld", value.positions.size());
    draw_label_opengl_widget(win, "norm", "%ld", value.normals.size());
    draw_label_opengl_widget(
        win, "texcoord", "%ld", value.texturecoords.size());
    draw_label_opengl_widget(win, "color", "%ld", value.colors.size());
    draw_label_opengl_widget(win, "radius", "%ld", value.radius.size());
    draw_label_opengl_widget(win, "tangsp", "%ld", value.tangentspaces.size());
    return edited;
}

inline bool draw_opengl_widgets_scene_inspector(
    const opengl_window& win, yocto_instance& value, yocto_scene& scene) {
    auto edited = 0;
    edited += draw_textinput_opengl_widget(win, "name", value.name);
    edited += draw_slider_opengl_widget(win, "frame[0]", value.frame.x, -1, 1);
    edited += draw_slider_opengl_widget(win, "frame[1]", value.frame.y, -1, 1);
    edited += draw_slider_opengl_widget(win, "frame[2]", value.frame.z, -1, 1);
    edited += draw_slider_opengl_widget(win, "frame.o", value.frame.o, -10, 10);
    edited += draw_combobox_opengl_widget(
        win, "shape", value.shape, scene.shapes, true);
    edited += draw_combobox_opengl_widget(
        win, "material", value.material, scene.materials, true);
    return edited;
}

inline bool draw_opengl_widgets_scene_inspector(
    const opengl_window& win, yocto_environment& value, yocto_scene& scene) {
    auto edited = 0;
    edited += draw_textinput_opengl_widget(win, "name", value.name);
    edited += draw_slider_opengl_widget(win, "frame[0]", value.frame.x, -1, 1);
    edited += draw_slider_opengl_widget(win, "frame[1]", value.frame.y, -1, 1);
    edited += draw_slider_opengl_widget(win, "frame[2]", value.frame.z, -1, 1);
    edited += draw_slider_opengl_widget(win, "frame.o", value.frame.o, -10, 10);
    edited += draw_hdr_coloredit_opengl_widget(win, "ke", value.emission);
    edited += draw_combobox_opengl_widget(
        win, "ke texture", value.emission_texture, scene.textures, true);
    return edited;
}

inline bool draw_opengl_widgets_scene_inspector(
    const opengl_window& win, yocto_scene_node& value, yocto_scene& scene) {
    auto edited = 0;
    edited += draw_textinput_opengl_widget(win, "name", value.name);
    edited += draw_combobox_opengl_widget(
        win, "parent", value.parent, scene.nodes, true);
    edited += draw_slider_opengl_widget(win, "local[0]", value.local.x, -1, 1);
    edited += draw_slider_opengl_widget(win, "local[1]", value.local.y, -1, 1);
    edited += draw_slider_opengl_widget(win, "local[2]", value.local.z, -1, 1);
    edited += draw_slider_opengl_widget(win, "local.o", value.local.o, -10, 10);
    edited += draw_slider_opengl_widget(
        win, "translation", value.translation, -10, 10);
    edited += draw_slider_opengl_widget(win, "rotation", value.rotation, -1, 1);
    edited += draw_slider_opengl_widget(win, "scale", value.scale, 0, 10);
    edited += draw_combobox_opengl_widget(
        win, "camera", value.camera, scene.cameras, true);
    edited += draw_combobox_opengl_widget(
        win, "instance", value.instance, scene.instances, true);
    edited += draw_combobox_opengl_widget(
        win, "environment", value.environment, scene.environments, true);
    return edited;
}

inline bool draw_opengl_widgets_scene_inspector(
    const opengl_window& win, yocto_animation& value, yocto_scene& scene) {
    auto edited = 0;
    edited += draw_textinput_opengl_widget(win, "name", value.name);
    edited += draw_textinput_opengl_widget(win, "path", value.filename);
    edited += draw_textinput_opengl_widget(win, "group", value.animation_group);
    // edited += draw_combobox_opengl_widget(win, "type", &value.type,
    // animation_type_names());
    draw_label_opengl_widget(win, "times", "%ld", value.keyframes_times.size());
    draw_label_opengl_widget(
        win, "translation", "%ld", value.translation_keyframes.size());
    draw_label_opengl_widget(
        win, "rotation", "%ld", value.rotation_keyframes.size());
    draw_label_opengl_widget(win, "scale", "%ld", value.scale_keyframes.size());
    draw_label_opengl_widget(
        win, "weights", "%ld", value.morph_weights_keyframes.size());
    draw_label_opengl_widget(win, "targets", "%ld", value.node_targets.size());
    return edited;
}

inline bool draw_opengl_widgets_scene_tree(const opengl_window& win,
    const string& lbl, yocto_scene& scene, pair<type_index, int>& sel,
    vector<pair<type_index, int>>& update_list, int height) {
    draw_opengl_widgets_scene_tree(win, scene, sel);
    auto update_len = update_list.size();
#if 0
    if (test_scn) {
        draw_add_elem_opengl_widgets(
            scene, "camera", scene.cameras, test_scn->cameras, sel, update_list);
        draw_add_elem_opengl_widgets(scene, "texture", scene.textures,
            test_scn->textures, sel, update_list);
        draw_add_elem_opengl_widgets(scene, "mat", scene.materials,
            test_scn->materials, sel, update_list);
        draw_add_elem_opengl_widgets(
            scene, "shape", scene.shapes, test_scn->shapes, sel, update_list);
        draw_add_elem_opengl_widgets(scene, "instance", scene.instances,
            test_scn->instances, sel, update_list);
        draw_add_elem_opengl_widgets(
            scene, "node", scene.nodes, test_scn->nodes, sel, update_list);
        draw_add_elem_opengl_widgets(scene, "environment", scene.environments,
            test_scn->environments, sel, update_list);
        draw_add_elem_opengl_widgets(scene, "anim", scene.animations,
            test_scn->animations, sel, update_list);
    }
#endif
    return update_list.size() != update_len;
}

inline bool draw_opengl_widgets_scene_inspector(const opengl_window& win,
    const string& lbl, yocto_scene& scene, pair<type_index, int>& sel,
    vector<pair<type_index, int>>& update_list, int height) {
    if (sel.first == typeid(void)) return false;

    auto update_len = update_list.size();

    if (sel.first == typeid(yocto_camera))
        if (draw_opengl_widgets_scene_inspector(
                win, scene.cameras[sel.second], scene))
            update_list.push_back({typeid(yocto_camera), sel.second});
    if (sel.first == typeid(yocto_shape))
        if (draw_opengl_widgets_scene_inspector(
                win, scene.shapes[sel.second], scene))
            update_list.push_back({typeid(yocto_shape), sel.second});
    if (sel.first == typeid(yocto_texture))
        if (draw_opengl_widgets_scene_inspector(
                win, scene.textures[sel.second], scene))
            update_list.push_back({typeid(yocto_texture), sel.second});
    if (sel.first == typeid(yocto_material))
        if (draw_opengl_widgets_scene_inspector(
                win, scene.materials[sel.second], scene))
            update_list.push_back({typeid(yocto_material), sel.second});
    if (sel.first == typeid(yocto_environment))
        if (draw_opengl_widgets_scene_inspector(
                win, scene.environments[sel.second], scene))
            update_list.push_back({typeid(yocto_environment), sel.second});
    if (sel.first == typeid(yocto_instance))
        if (draw_opengl_widgets_scene_inspector(
                win, scene.instances[sel.second], scene))
            update_list.push_back({typeid(yocto_instance), sel.second});
    if (sel.first == typeid(yocto_scene_node))
        if (draw_opengl_widgets_scene_inspector(
                win, scene.nodes[sel.second], scene))
            update_list.push_back({typeid(yocto_scene_node), sel.second});
    if (sel.first == typeid(yocto_animation))
        if (draw_opengl_widgets_scene_inspector(
                win, scene.animations[sel.second], scene))
            update_list.push_back({typeid(yocto_animation), sel.second});

    return update_list.size() != update_len;
}

#endif
