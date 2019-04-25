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

#include <any>
#include <typeindex>
#include <unordered_map>
using std::any;
using std::any_cast;
using std::type_index;
using std::unordered_map;

struct app_selection {
    type_index type  = typeid(void);
    int        index = -1;
};

struct app_edit {
    type_index type   = typeid(void);
    int        index  = -1;
    any        data   = {};
    bool       reload = false;
};

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
    app_selection& sel);

template <typename T>
inline void draw_opengl_widgets_scene_tree(const opengl_window& win,
    const string& lbl_, yocto_scene& scene, int index, const vector<T*>& vals,
    app_selection& sel);

inline void draw_scene_tree_opengl_widgets_rec(const opengl_window& win,
    const string& lbl_, yocto_scene& scene, const yocto_camera& value,
    app_selection& sel) {}
inline void draw_scene_tree_opengl_widgets_rec(const opengl_window& win,
    const string& lbl_, yocto_scene& scene, const yocto_texture& value,
    app_selection& sel) {}
inline void draw_scene_tree_opengl_widgets_rec(const opengl_window& win,
    const string& lbl_, yocto_scene& scene, const yocto_voltexture& value,
    app_selection& sel) {}
inline void draw_scene_tree_opengl_widgets_rec(const opengl_window& win,
    const string& lbl_, yocto_scene& scene, const yocto_material& value,
    app_selection& sel) {
    draw_opengl_widgets_scene_tree(
        win, "emission", scene, value.emission_texture, scene.textures, sel);
    draw_opengl_widgets_scene_tree(
        win, "diffuse", scene, value.diffuse_texture, scene.textures, sel);
    draw_opengl_widgets_scene_tree(
        win, "specular", scene, value.specular_texture, scene.textures, sel);
    draw_opengl_widgets_scene_tree(
        win, "normal", scene, value.normal_texture, scene.textures, sel);
}
inline void draw_scene_tree_opengl_widgets_rec(const opengl_window& win,
    const string& lbl_, yocto_scene& scene, const yocto_shape& value,
    app_selection& sel) {}
inline void draw_scene_tree_opengl_widgets_rec(const opengl_window& win,
    const string& lbl_, yocto_scene& scene, const yocto_subdiv& value,
    app_selection& sel) {
    draw_opengl_widgets_scene_tree(
        win, "shapes", scene, value.tesselated_shape, scene.shapes, sel);
    draw_opengl_widgets_scene_tree(win, "displament", scene,
        value.displacement_texture, scene.textures, sel);
}

inline void draw_scene_tree_opengl_widgets_rec(const opengl_window& win,
    const string& lbl_, yocto_scene& scene, const yocto_instance& value,
    app_selection& sel) {
    draw_opengl_widgets_scene_tree(
        win, "shape", scene, value.shape, scene.shapes, sel);
    draw_opengl_widgets_scene_tree(
        win, "material", scene, value.material, scene.materials, sel);
}
inline void draw_scene_tree_opengl_widgets_rec(const opengl_window& win,
    const string& lbl_, yocto_scene& scene, const yocto_environment& value,
    app_selection& sel) {
    draw_opengl_widgets_scene_tree(
        win, "emission", scene, value.emission_texture, scene.textures, sel);
}
inline void draw_scene_tree_opengl_widgets_rec(const opengl_window& win,
    const string& lbl_, yocto_scene& scene, const yocto_scene_node& value,
    app_selection& sel) {
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
    app_selection& sel) {
    auto tid = 0;
    for (auto tg : value.node_targets) {
        draw_opengl_widgets_scene_tree(
            win, "target" + std::to_string(tid++), scene, tg, scene.nodes, sel);
    }
}

template <typename T>
inline void draw_opengl_widgets_scene_tree(const opengl_window& win,
    const string& lbl_, yocto_scene& scene, int index, const vector<T>& vals,
    app_selection& sel) {
    if (index < 0) return;
    auto lbl = vals[index].uri;
    if (!empty(lbl_)) lbl = lbl_ + ": " + vals[index].uri;
    auto selected = sel.type == type_index(typeid(T)) && sel.index == index;
    if (begin_selectabletreenode_opengl_widget(win, lbl.c_str(), selected)) {
        draw_scene_tree_opengl_widgets_rec(win, lbl_, scene, vals[index], sel);
        end_treenode_opengl_widget(win);
    }
    if (selected) sel = {type_index(typeid(T)), index};
}

template <typename T>
inline void draw_opengl_widgets_scene_tree(const opengl_window& win,
    const string& lbl_, yocto_scene& scene, int index, const vector<T*>& vals,
    app_selection& sel) {
    if (index < 0) return;
    auto lbl = vals[index]->name;
    if (!empty(lbl_)) lbl = lbl_ + ": " + vals[index]->name;
    auto selected = sel.type == type_index(typeid(T)) && sel.index == index;
    if (begin_selectabletreenode_opengl_widget(win, lbl.c_str(), selected)) {
        draw_scene_tree_opengl_widgets_rec(win, lbl_, scene, vals[index], sel);
        end_treenode_opengl_widget(win);
    }
    if (selected) sel = {type_index(typeid(T)), index};
}

inline void draw_opengl_widgets_scene_tree(
    const opengl_window& win, yocto_scene& scene, app_selection& sel) {
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
inline bool draw_opengl_widgets_scene_inspector(const opengl_window& win,
    const yocto_camera& value, const app_selection& sel, app_edit& edit,
    yocto_scene& scene) {
    auto edited_value = value;
    auto edited       = false;
    if (draw_textinput_opengl_widget(win, "uri", edited_value.uri)) {
        edited = true;
    }
    if (draw_slider_opengl_widget(
            win, "frame.x", edited_value.frame.x, -1, 1)) {
        edited = true;
    }
    if (draw_slider_opengl_widget(
            win, "frame.y", edited_value.frame.y, -1, 1)) {
        edited = true;
    }
    if (draw_slider_opengl_widget(
            win, "frame.z", edited_value.frame.z, -1, 1)) {
        edited = true;
    }
    if (draw_slider_opengl_widget(
            win, "frame.o", edited_value.frame.o, -10, 10)) {
        edited = true;
    }
    if (draw_checkbox_opengl_widget(win, "ortho", edited_value.orthographic)) {
        edited = true;
    }
    if (draw_slider_opengl_widget(
            win, "film width", edited_value.film_width, 0.01f, 1)) {
        edited = true;
    }
    if (draw_slider_opengl_widget(
            win, "film height", edited_value.film_height, 0.01f, 1)) {
        edited = true;
    }
    if (draw_slider_opengl_widget(
            win, "focal length", edited_value.focal_length, 0.01f, 1)) {
        edited = true;
    }
    if (draw_slider_opengl_widget(
            win, "focus distance", edited_value.focus_distance, 0.01f, 1000)) {
        edited = true;
    }
    if (draw_slider_opengl_widget(
            win, "lens aperture", edited_value.lens_aperture, 0, 5)) {
        edited = true;
    }
    auto from = edited_value.frame.o,
         to   = edited_value.frame.o -
              edited_value.focus_distance * edited_value.frame.z;
    draw_slider_opengl_widget(win, "!!from", from, -10, 10);
    draw_slider_opengl_widget(win, "!!to", to, -10, 10);
    if (edited) {
        edit = {sel.type, sel.index, edited_value, false};
    }
    return edited;
}

/// Visit struct elements.
inline bool draw_opengl_widgets_scene_inspector(const opengl_window& win,
    const yocto_texture& value, const app_selection& sel, app_edit& edit,
    yocto_scene& scene) {
    auto edited_value = yocto_texture{};
    edited_value.uri  = value.uri;
    auto edited       = false;
    if (draw_textinput_opengl_widget(win, "uri", edited_value.uri)) {
        edited = true;
    }
    draw_label_opengl_widget(win, "hdr_image", "%d x %d",
        value.hdr_image.size().x, value.hdr_image.size().y);
    draw_label_opengl_widget(win, "ldr_image", "%d x %d",
        value.ldr_image.size().x, value.ldr_image.size().y);
    if (edited) {
        auto reload = edited_value.uri != value.uri;
        if (!reload) {
            edited_value.hdr_image = value.hdr_image;
            edited_value.ldr_image = value.ldr_image;
        }
        edit = {sel.type, sel.index, edited_value, reload};
    }
    return edited;
}

inline bool draw_opengl_widgets_scene_inspector(const opengl_window& win,
    const yocto_voltexture& value, const app_selection& sel, app_edit& edit,
    yocto_scene& scene) {
    auto edited_value = yocto_voltexture{};
    edited_value.uri  = value.uri;
    auto edited       = false;
    if (draw_textinput_opengl_widget(win, "uri", edited_value.uri)) {
        edited = true;
    }
    draw_label_opengl_widget(win, "voxel_data", "%d x %d x %d",
        value.volume_data.size().x, value.volume_data.size().y,
        value.volume_data.size().z);
    if (edited) {
        auto reload = edited_value.uri != value.uri;
        if (!reload) {
            edited_value.volume_data = value.volume_data;
        }
        edit = {sel.type, sel.index, edited_value, reload};
    }
    return edited;
}

inline bool draw_opengl_widgets_scene_inspector(const opengl_window& win,
    const yocto_material& value, const app_selection& sel, app_edit& edit,
    yocto_scene& scene) {
    auto edited_value = value;
    auto edited       = false;
    if (draw_textinput_opengl_widget(win, "uri", edited_value.uri)) {
        edited = true;
    }
    if (draw_hdr_coloredit_opengl_widget(
            win, "emission", edited_value.emission)) {
        edited = true;
    }
    if (draw_coloredit_opengl_widget(win, "diffuse", edited_value.diffuse)) {
        edited = true;
    }
    if (draw_coloredit_opengl_widget(win, "specular", edited_value.specular)) {
        edited = true;
    }
    if (draw_coloredit_opengl_widget(
            win, "transmission", edited_value.transmission)) {
        edited = true;
    }
    if (draw_slider_opengl_widget(
            win, "roughness", edited_value.roughness, 0, 1)) {
        edited = true;
    }
    if (draw_slider_opengl_widget(win, "opacity", edited_value.opacity, 0, 1)) {
        edited = true;
    }
    if (draw_checkbox_opengl_widget(win, "fresnel", edited_value.fresnel)) {
        edited = true;
    }
    continue_opengl_widget_line(win);
    if (draw_checkbox_opengl_widget(win, "refract", edited_value.refract)) {
        edited = true;
    }
    if (draw_coloredit_opengl_widget(
            win, "volume_density", edited_value.volume_density)) {
        edited = true;
    }  // 0, 10
    if (draw_coloredit_opengl_widget(
            win, "volume_albedo", edited_value.volume_albedo)) {
        edited = true;
    }  // 0, 1
    if (draw_slider_opengl_widget(
            win, "volume_phaseg", edited_value.volume_phaseg, -1, 1)) {
        edited = true;
    }
    if (draw_combobox_opengl_widget(win, "emission_texture",
            edited_value.emission_texture, scene.textures, true)) {
        edited = true;
    }
    if (draw_combobox_opengl_widget(win, "diffuse_texture",
            edited_value.diffuse_texture, scene.textures, true)) {
        edited = true;
    }
    if (draw_combobox_opengl_widget(win, "specular_texture",
            edited_value.specular_texture, scene.textures, true)) {
        edited = true;
    }
    if (draw_combobox_opengl_widget(win, "transmission_texture",
            edited_value.transmission_texture, scene.textures, true)) {
        edited = true;
    }
    if (draw_combobox_opengl_widget(win, "roughness_texture",
            edited_value.roughness_texture, scene.textures, true)) {
        edited = true;
    }
    if (draw_combobox_opengl_widget(win, "normal_texture",
            edited_value.normal_texture, scene.textures, true)) {
        edited = true;
    }
    if (draw_checkbox_opengl_widget(
            win, "base metallic", edited_value.base_metallic)) {
        edited = true;
    }
    if (draw_checkbox_opengl_widget(
            win, "glTF textures", edited_value.gltf_textures)) {
        edited = true;
    }
    if (edited) {
        edit = {sel.type, sel.index, edited_value, false};
    }
    return edited;
}

inline bool draw_opengl_widgets_scene_inspector(const opengl_window& win,
    const yocto_shape& value, const app_selection& sel, app_edit& edit,
    yocto_scene& scene) {
    auto edited_value = yocto_shape{};
    edited_value.uri  = value.uri;
    auto edited       = false;
    if (draw_textinput_opengl_widget(win, "uri", edited_value.uri)) {
        edited = true;
    }
    draw_label_opengl_widget(win, "points", "%ld", value.points.size());
    draw_label_opengl_widget(win, "lines", "%ld", value.lines.size());
    draw_label_opengl_widget(win, "triangles", "%ld", value.triangles.size());
    draw_label_opengl_widget(win, "quads", "%ld", value.quads.size());
    draw_label_opengl_widget(
        win, "quads pos", "%ld", value.quads_positions.size());
    draw_label_opengl_widget(
        win, "quads norm", "%ld", value.quads_normals.size());
    draw_label_opengl_widget(
        win, "quads texcoord", "%ld", value.quads_texcoords.size());
    draw_label_opengl_widget(win, "pos", "%ld", value.positions.size());
    draw_label_opengl_widget(win, "norm", "%ld", value.normals.size());
    draw_label_opengl_widget(win, "texcoord", "%ld", value.texcoords.size());
    draw_label_opengl_widget(win, "color", "%ld", value.colors.size());
    draw_label_opengl_widget(win, "radius", "%ld", value.radius.size());
    draw_label_opengl_widget(win, "tangsp", "%ld", value.tangents.size());
    if (edited) {
        auto reload = edited_value.uri != value.uri;
        if (!reload) {
            edited_value.points          = value.points;
            edited_value.lines           = value.lines;
            edited_value.triangles       = value.triangles;
            edited_value.quads           = value.quads;
            edited_value.quads_positions = value.quads_positions;
            edited_value.quads_normals   = value.quads_normals;
            edited_value.quads_texcoords = value.quads_texcoords;
            edited_value.positions       = value.positions;
            edited_value.normals         = value.normals;
            edited_value.texcoords       = value.texcoords;
            edited_value.colors          = value.colors;
            edited_value.radius          = value.radius;
            edited_value.tangents        = value.tangents;
        }
        edit = {sel.type, sel.index, edited_value, reload};
    }
    return edited;
    return edited;
}

inline bool draw_opengl_widgets_scene_inspector(const opengl_window& win,
    const yocto_subdiv& value, const app_selection& sel, app_edit& edit,
    yocto_scene& scene) {
    auto edited_value                 = yocto_subdiv{};
    edited_value.uri                  = value.uri;
    edited_value.subdivision_level    = value.subdivision_level;
    edited_value.catmull_clark        = value.catmull_clark;
    edited_value.compute_normals      = value.compute_normals;
    edited_value.preserve_facevarying = value.preserve_facevarying;
    edited_value.tesselated_shape     = value.tesselated_shape;
    edited_value.displacement_texture = value.displacement_texture;
    edited_value.displacement_scale   = value.displacement_scale;
    auto edited                       = false;
    if (draw_textinput_opengl_widget(win, "uri", edited_value.uri)) {
        edited = true;
    }
    if (draw_slider_opengl_widget(
            win, "subdivision_level", edited_value.subdivision_level, 0, 10)) {
        edited = true;
    }
    if (draw_checkbox_opengl_widget(
            win, "catmull_clark", edited_value.catmull_clark)) {
        edited = true;
    }
    if (draw_checkbox_opengl_widget(
            win, "compute_normals", edited_value.compute_normals)) {
        edited = true;
    }
    if (draw_checkbox_opengl_widget(
            win, "preserve_facevarying", edited_value.preserve_facevarying)) {
        edited = true;
    }
    if (draw_combobox_opengl_widget(win, "tesselated_shape",
            edited_value.tesselated_shape, scene.textures, true)) {
        edited = true;
    }
    if (draw_combobox_opengl_widget(win, "displacement_texture",
            edited_value.displacement_texture, scene.textures, true)) {
        edited = true;
    }
    if (draw_slider_opengl_widget(
            win, "displacement_scale", edited_value.displacement_scale, 0, 1)) {
        edited = true;
    }
    draw_label_opengl_widget(win, "points", "%ld", value.points.size());
    draw_label_opengl_widget(win, "lines", "%ld", value.lines.size());
    draw_label_opengl_widget(win, "triangles", "%ld", value.triangles.size());
    draw_label_opengl_widget(win, "quads", "%ld", value.quads.size());
    draw_label_opengl_widget(
        win, "quads pos", "%ld", value.quads_positions.size());
    draw_label_opengl_widget(
        win, "quads norm", "%ld", value.quads_normals.size());
    draw_label_opengl_widget(
        win, "quads texcoord", "%ld", value.quads_texcoords.size());
    draw_label_opengl_widget(win, "pos", "%ld", value.positions.size());
    draw_label_opengl_widget(win, "norm", "%ld", value.normals.size());
    draw_label_opengl_widget(win, "texcoord", "%ld", value.texcoords.size());
    draw_label_opengl_widget(win, "color", "%ld", value.colors.size());
    draw_label_opengl_widget(win, "radius", "%ld", value.radius.size());
    if (edited) {
        auto reload = edited_value.uri != value.uri;
        if (!reload) {
            edited_value.points          = value.points;
            edited_value.lines           = value.lines;
            edited_value.triangles       = value.triangles;
            edited_value.quads           = value.quads;
            edited_value.quads_positions = value.quads_positions;
            edited_value.quads_normals   = value.quads_normals;
            edited_value.quads_texcoords = value.quads_texcoords;
            edited_value.positions       = value.positions;
            edited_value.normals         = value.normals;
            edited_value.texcoords       = value.texcoords;
            edited_value.colors          = value.colors;
            edited_value.radius          = value.radius;
        }
        edit = {sel.type, sel.index, edited_value, reload};
    }
    return edited;
}

inline bool draw_opengl_widgets_scene_inspector(const opengl_window& win,
    const yocto_instance& value, const app_selection& sel, app_edit& edit,
    yocto_scene& scene) {
    auto edited_value = value;
    auto edited       = false;
    if (draw_textinput_opengl_widget(win, "uri", edited_value.uri)) {
        edited = true;
    }
    if (draw_slider_opengl_widget(
            win, "frame[0]", edited_value.frame.x, -1, 1)) {
        edited = true;
    }
    if (draw_slider_opengl_widget(
            win, "frame[1]", edited_value.frame.y, -1, 1)) {
        edited = true;
    }
    if (draw_slider_opengl_widget(
            win, "frame[2]", edited_value.frame.z, -1, 1)) {
        edited = true;
    }
    if (draw_slider_opengl_widget(
            win, "frame.o", edited_value.frame.o, -10, 10)) {
        edited = true;
    }
    if (draw_combobox_opengl_widget(
            win, "shape", edited_value.shape, scene.shapes, true)) {
        edited = true;
    }
    if (draw_combobox_opengl_widget(
            win, "material", edited_value.material, scene.materials, true)) {
        edited = true;
    }
    if (edited) {
        edit = {sel.type, sel.index, edited_value, false};
    }
    return edited;
}

inline bool draw_opengl_widgets_scene_inspector(const opengl_window& win,
    const yocto_environment& value, const app_selection& sel, app_edit& edit,
    yocto_scene& scene) {
    auto edited_value = value;
    auto edited       = false;
    if (draw_textinput_opengl_widget(win, "uri", edited_value.uri)) {
        edited = true;
    }
    if (draw_slider_opengl_widget(
            win, "frame[0]", edited_value.frame.x, -1, 1)) {
        edited = true;
    }
    if (draw_slider_opengl_widget(
            win, "frame[1]", edited_value.frame.y, -1, 1)) {
        edited = true;
    }
    if (draw_slider_opengl_widget(
            win, "frame[2]", edited_value.frame.z, -1, 1)) {
        edited = true;
    }
    if (draw_slider_opengl_widget(
            win, "frame.o", edited_value.frame.o, -10, 10)) {
        edited = true;
    }
    if (draw_hdr_coloredit_opengl_widget(win, "ke", edited_value.emission)) {
        edited = true;
    }
    if (draw_combobox_opengl_widget(win, "ke texture",
            edited_value.emission_texture, scene.textures, true)) {
        edited = true;
    }
    if (edited) {
        edit = {sel.type, sel.index, edited_value, false};
    }
    return edited;
}

inline bool draw_opengl_widgets_scene_inspector(const opengl_window& win,
    const yocto_scene_node& value, const app_selection& sel, app_edit& edit,
    yocto_scene& scene) {
    auto edited_value = value;
    auto edited       = false;
    if (draw_textinput_opengl_widget(win, "uri", edited_value.uri)) {
        edited = true;
    }
    if (draw_combobox_opengl_widget(
            win, "parent", edited_value.parent, scene.nodes, true)) {
        edited = true;
    }
    if (draw_slider_opengl_widget(
            win, "local[0]", edited_value.local.x, -1, 1)) {
        edited = true;
    }
    if (draw_slider_opengl_widget(
            win, "local[1]", edited_value.local.y, -1, 1)) {
        edited = true;
    }
    if (draw_slider_opengl_widget(
            win, "local[2]", edited_value.local.z, -1, 1)) {
        edited = true;
    }
    if (draw_slider_opengl_widget(
            win, "local.o", edited_value.local.o, -10, 10)) {
        edited = true;
    }
    if (draw_slider_opengl_widget(
            win, "translation", edited_value.translation, -10, 10)) {
        edited = true;
    }
    if (draw_slider_opengl_widget(
            win, "rotation", edited_value.rotation, -1, 1)) {
        edited = true;
    }
    if (draw_slider_opengl_widget(win, "scale", edited_value.scale, 0, 10)) {
        edited = true;
    }
    if (draw_combobox_opengl_widget(
            win, "camera", edited_value.camera, scene.cameras, true)) {
        edited = true;
    }
    if (draw_combobox_opengl_widget(
            win, "instance", edited_value.instance, scene.instances, true)) {
        edited = true;
    }
    if (draw_combobox_opengl_widget(win, "environment",
            edited_value.environment, scene.environments, true)) {
        edited = true;
    }
    if (edited) {
        edit = {sel.type, sel.index, edited_value, false};
    }
    return edited;
}

inline bool draw_opengl_widgets_scene_inspector(const opengl_window& win,
    const yocto_animation& value, const app_selection& sel, app_edit& edit,
    yocto_scene& scene) {
    auto edited_value = value;
    auto edited       = false;
    if (draw_textinput_opengl_widget(win, "uri", edited_value.uri)) {
        edited = true;
    }
    if (draw_textinput_opengl_widget(win, "path", edited_value.filename)) {
        edited = true;
    }
    if (draw_textinput_opengl_widget(
            win, "group", edited_value.animation_group)) {
        edited = true;
    }
    // if(draw_combobox_opengl_widget(win, "type", &value.type,
    // animation_type_names())) edited = false;
    draw_label_opengl_widget(
        win, "times", "%ld", edited_value.keyframes_times.size());
    draw_label_opengl_widget(
        win, "translation", "%ld", edited_value.translation_keyframes.size());
    draw_label_opengl_widget(
        win, "rotation", "%ld", edited_value.rotation_keyframes.size());
    draw_label_opengl_widget(
        win, "scale", "%ld", edited_value.scale_keyframes.size());
    draw_label_opengl_widget(
        win, "weights", "%ld", value.morph_weights_keyframes.size());
    draw_label_opengl_widget(
        win, "targets", "%ld", edited_value.node_targets.size());
    if (edited) {
        edit = {sel.type, sel.index, edited_value, false};
    }
    return edited;
}

inline void draw_opengl_widgets_scene_tree(const opengl_window& win,
    const string& lbl, yocto_scene& scene, app_selection& sel, int height) {
    draw_opengl_widgets_scene_tree(win, scene, sel);
}

inline bool draw_opengl_widgets_scene_inspector(const opengl_window& win,
    const string& lbl, yocto_scene& scene, app_selection& sel, app_edit& edit,
    int height) {
    if (sel.type == typeid(void)) return false;

    edit = app_edit{};
    if (sel.type == typeid(yocto_camera)) {
        return draw_opengl_widgets_scene_inspector(
            win, scene.cameras[sel.index], sel, edit, scene);
    }
    if (sel.type == typeid(yocto_shape)) {
        return draw_opengl_widgets_scene_inspector(
            win, scene.shapes[sel.index], sel, edit, scene);
    }
    if (sel.type == typeid(yocto_subdiv)) {
        return draw_opengl_widgets_scene_inspector(
            win, scene.subdivs[sel.index], sel, edit, scene);
    }
    if (sel.type == typeid(yocto_texture)) {
        return draw_opengl_widgets_scene_inspector(
            win, scene.textures[sel.index], sel, edit, scene);
    }
    if (sel.type == typeid(yocto_voltexture)) {
        return draw_opengl_widgets_scene_inspector(
            win, scene.voltextures[sel.index], sel, edit, scene);
    }
    if (sel.type == typeid(yocto_material)) {
        return draw_opengl_widgets_scene_inspector(
            win, scene.materials[sel.index], sel, edit, scene);
    }
    if (sel.type == typeid(yocto_environment)) {
        return draw_opengl_widgets_scene_inspector(
            win, scene.environments[sel.index], sel, edit, scene);
    }
    if (sel.type == typeid(yocto_instance)) {
        return draw_opengl_widgets_scene_inspector(
            win, scene.instances[sel.index], sel, edit, scene);
    }
    if (sel.type == typeid(yocto_scene_node)) {
        return draw_opengl_widgets_scene_inspector(
            win, scene.nodes[sel.index], sel, edit, scene);
    }
    if (sel.type == typeid(yocto_animation)) {
        return draw_opengl_widgets_scene_inspector(
            win, scene.animations[sel.index], sel, edit, scene);
    }

    return false;
}

#endif
