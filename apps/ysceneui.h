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
inline void draw_scene_tree_glwidgets_rec(const std::string& lbl_, T* val,
    void*& sel) {}

template <typename T>
inline void draw_glwidgets_scene_tree(const std::string& lbl_, T* val,
    void*& sel) {
    if (!val) return;
    auto lbl = val->name;
    if (!lbl_.empty()) lbl = lbl_ + ": " + val->name;
    if(ImGui::SelectableTreeNode(lbl.c_str(), &sel, val)) {
        draw_scene_tree_glwidgets_rec(lbl_, val, sel);
        ImGui::TreePop();
    }
}

template <>
inline void draw_scene_tree_glwidgets_rec<ygl::instance>(const std::string& lbl_,
    ygl::instance* val, void*& sel) {
    draw_glwidgets_scene_tree("shp", val->shp, sel);
    draw_glwidgets_scene_tree("sbd", val->sbd, sel);
    draw_glwidgets_scene_tree("mat", val->mat, sel);
}

template <>
inline void draw_scene_tree_glwidgets_rec<ygl::material>(const std::string& lbl_,
    ygl::material* val, void*& sel) {
    draw_glwidgets_scene_tree("ke", val->ke_txt, sel);
    draw_glwidgets_scene_tree("kd", val->kd_txt, sel);
    draw_glwidgets_scene_tree("ks", val->ks_txt, sel);
    draw_glwidgets_scene_tree("bump", val->bump_txt, sel);
    draw_glwidgets_scene_tree("disp", val->disp_txt, sel);
    draw_glwidgets_scene_tree("norm", val->norm_txt, sel);
}
template <>
inline void draw_scene_tree_glwidgets_rec<ygl::environment>(const std::string& lbl_,
    ygl::environment* val, void*& sel) {
    draw_glwidgets_scene_tree("ke", val->ke_txt, sel);
}
template <>
inline void draw_scene_tree_glwidgets_rec<ygl::node>(const std::string& lbl_,
    ygl::node* val, void*& sel) {
    draw_glwidgets_scene_tree("ist", val->ist, sel);
    draw_glwidgets_scene_tree("cam", val->cam, sel);
    draw_glwidgets_scene_tree("env", val->env, sel);
    draw_glwidgets_scene_tree("par", val->parent, sel);
    auto cid = 0;
    for (auto ch : val->children) {
        draw_glwidgets_scene_tree(
            "ch" + std::to_string(cid++), ch, sel);
    }
}
template <>
inline void draw_scene_tree_glwidgets_rec<ygl::animation>(const std::string& lbl_,
    ygl::animation* val, void*& sel) {
    auto tid = 0;
    for (auto tg : val->targets) {
        draw_glwidgets_scene_tree(
            "tg" + std::to_string(tid++), tg, sel);
    }
}

inline void draw_glwidgets_scene_tree(ygl::scene* scn, void*& sel) {
    if (!scn->cameras.empty() && ImGui::TreeNode("cameras")) {
        for (auto v : scn->cameras)
            draw_glwidgets_scene_tree("", v, sel);
        ImGui::TreePop();
    }
    if (!scn->shapes.empty() && ImGui::TreeNode("shapes")) {
        for (auto v : scn->shapes)
            draw_glwidgets_scene_tree("", v, sel);
        ImGui::TreePop();
    }
    if (!scn->subdivs.empty() && ImGui::TreeNode("subdivs")) {
        for (auto v : scn->subdivs)
            draw_glwidgets_scene_tree("", v, sel);
        ImGui::TreePop();
    }
    if (!scn->instances.empty() && ImGui::TreeNode("instances")) {
        for (auto v : scn->instances)
            draw_glwidgets_scene_tree("", v, sel);
        ImGui::TreePop();
    }
    if (!scn->materials.empty() && ImGui::TreeNode("materials")) {
        for (auto v : scn->materials)
            draw_glwidgets_scene_tree("", v, sel);
        ImGui::TreePop();
    }
    if (!scn->textures.empty() && ImGui::TreeNode("textures")) {
        for (auto v : scn->textures)
            draw_glwidgets_scene_tree("", v, sel);
        ImGui::TreePop();
    }
    if (!scn->environments.empty() && ImGui::TreeNode("environments")) {
        for (auto v : scn->environments)
            draw_glwidgets_scene_tree("", v, sel);
        ImGui::TreePop();
    }
    if (!scn->nodes.empty() && ImGui::TreeNode("nodes")) {
        for (auto v : scn->nodes)
            draw_glwidgets_scene_tree("", v, sel);
        ImGui::TreePop();
    }
    if (!scn->animations.empty() && ImGui::TreeNode("animations")) {
        for (auto v : scn->animations)
            draw_glwidgets_scene_tree("", v, sel);
        ImGui::TreePop();
    }
}

/// Visit struct elements.
inline bool draw_glwidgets_scene_inspector(ygl::camera* val, ygl::scene* scn) {
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
inline bool draw_glwidgets_scene_inspector(ygl::texture* val, ygl::scene* scn) {
    auto edited = 0;
    edited += ImGui::InputText("name", &val->name);
    edited += ImGui::InputText("path", &val->path);
    edited += ImGui::Checkbox("clamp", &val->clamp);
    edited += ImGui::SliderFloat("scale", &val->scale, 0, 1);
    edited += ImGui::SliderFloat("gamma", &val->gamma, 1, 2.2f);
    ImGui::LabelText("img", "%d x %d", val->img.width, val->img.height);
    return edited;
}

inline bool draw_glwidgets_scene_inspector(ygl::material* val, ygl::scene* scn) {
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

inline bool draw_glwidgets_scene_inspector(ygl::shape* val, ygl::scene* scn) {
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

inline bool draw_glwidgets_scene_inspector(ygl::subdiv* val, ygl::scene* scn) {
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

inline bool draw_glwidgets_scene_inspector(ygl::instance* val, ygl::scene* scn) {
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

inline bool draw_glwidgets_scene_inspector(ygl::environment* val, ygl::scene* scn) {
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

inline bool draw_glwidgets_scene_inspector(ygl::node* val, ygl::scene* scn) {
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

inline bool draw_glwidgets_scene_inspector(ygl::animation* val, ygl::scene* scn) {
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

inline bool draw_glwidgets_scene_tree(const std::string& lbl, ygl::scene* scn,
    void*& sel, std::vector<std::pair<std::string, void*>>& update_list,
    int height) {
    if (!scn) return false;
    ImGui::PushID(scn);
    ImGui::BeginChild("scrolling scene tree", ImVec2(0, height), false);
    draw_glwidgets_scene_tree(scn, sel);

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

inline bool draw_glwidgets_scene_inspector(const std::string& lbl, ygl::scene* scn,
    void*& sel, std::vector<std::pair<std::string, void*>>& update_list,
    int height) {
    if (!scn || !sel) return false;
    ImGui::PushID(sel);
    ImGui::BeginChild("scrolling scene inspector", ImVec2(0, height), false);

    auto update_len = update_list.size();

    for(auto cam : scn->cameras) {
        if(cam != sel) continue;
        if(draw_glwidgets_scene_inspector(cam, scn))
            update_list.push_back({"camera", cam});
    }
    for(auto shp : scn->shapes) {
        if(shp != sel) continue ;
        if(draw_glwidgets_scene_inspector(shp, scn))
            update_list.push_back({"shape", shp});
    }
    for(auto sbd : scn->subdivs) {
        if(sbd != sel) continue ;
        if(draw_glwidgets_scene_inspector(sbd, scn))
            update_list.push_back({"subdiv", sbd});
    }
    for(auto txt : scn->textures) {
        if(txt != sel) continue ;
        if(draw_glwidgets_scene_inspector(txt, scn))
            update_list.push_back({"texture", txt});
    }
    for(auto mat : scn->materials) {
        if(mat != sel) continue ;
        if(draw_glwidgets_scene_inspector(mat, scn))
            update_list.push_back({"material", mat});
    }
    for(auto env : scn->environments) {
        if(env != sel) continue ;
        if(draw_glwidgets_scene_inspector(env, scn))
            update_list.push_back({"environment", env});
    }
    for(auto ist : scn->instances) {
        if(ist != sel) continue ;
        if(draw_glwidgets_scene_inspector(ist, scn))
            update_list.push_back({"instance", ist});
    }
    for(auto nde : scn->nodes) {
        if(nde != sel) continue ;
        if(draw_glwidgets_scene_inspector(nde, scn))
            update_list.push_back({"node", nde});
    }
    for(auto anm : scn->animations) {
        if(anm != sel) continue ;
        if(draw_glwidgets_scene_inspector(anm, scn))
            update_list.push_back({"animation", anm});
    }

    ImGui::EndChild();
    ImGui::PopID();
    return update_list.size() != update_len;
}

#endif
