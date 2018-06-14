//
// Implementation of yapp_ui.h
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2018 Fabio Pellacini
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//

#include "yapp_ui.h"

namespace ygl {

static const std::map<animation_type, std::string>& animation_type_names() {
    static auto names = std::map<animation_type, std::string>{
        {animation_type::linear, "linear"},
        {animation_type::step, "step"},
        {animation_type::bezier, "bezier"},
    };
    return names;
}

// Vertex buffers for scene drawing. Members are not part of the public API.
struct glshape {
    glvertex_buffer pos;        // position
    glvertex_buffer norm;       // normals
    glvertex_buffer texcoord;   // texcoord
    glvertex_buffer texcoord1;  // texcoord
    glvertex_buffer tangsp;     // tangent space
    glvertex_buffer color;      // color
    glvertex_buffer points;     // point elements
    glvertex_buffer lines;      // line elements
    glvertex_buffer triangles;  // triangle elements
    glvertex_buffer quads;      // quad elements as tris
    glvertex_buffer beziers;    // bezier elements as l.
    glvertex_buffer edges;      // edge elements
};

void update_gldata(const std::shared_ptr<texture>& txt) {
    if (!txt->gl_data) txt->gl_data = new gltexture();
    auto& gtxt = *(gltexture*)txt->gl_data;
    if (!txt->img.empty()) {
        update_gltexture(gtxt, txt->img.width(), txt->img.height(),
            txt->img.pixels(), true, true, false, false);
    }
}
void update_gldata(const std::shared_ptr<shape>& shp) {
    auto update_vert_buffer = [](auto& buf, const auto& vert) {
        if (vert.empty()) {
            clear_glbuffer(buf);
        } else {
            update_glbuffer(buf, false, vert);
        }
    };
    auto update_elem_buffer = [](auto& buf, const auto& elem) {
        if (elem.empty()) {
            clear_glbuffer(buf);
        } else {
            update_glbuffer(buf, true, elem);
        }
    };

    if (!shp->gl_data) shp->gl_data = new glshape();
    auto& gshp = *(glshape*)shp->gl_data;
    update_vert_buffer(gshp.pos, shp->pos);
    update_vert_buffer(gshp.norm, shp->norm);
    update_vert_buffer(gshp.texcoord, shp->texcoord);
    update_vert_buffer(gshp.color, shp->color);
    update_vert_buffer(gshp.tangsp, shp->tangsp);
    update_elem_buffer(gshp.points, shp->points);
    update_elem_buffer(gshp.lines, shp->lines);
    update_elem_buffer(gshp.triangles, shp->triangles);
    update_elem_buffer(gshp.edges, get_edges(shp->triangles));
}
void update_gldata(const std::shared_ptr<scene>& scn) {
    for (auto txt : scn->textures) update_gldata(txt);
    for (auto shp : scn->shapes) update_gldata(shp);
}
void clear_gldata(const std::shared_ptr<texture>& txt) {
    if (!txt->gl_data) return;
    clear_gltexture(*(gltexture*)txt->gl_data);
    delete (gltexture*)txt->gl_data;
}
void clear_gldata(const std::shared_ptr<shape>& shp) {
    if (!shp->gl_data) return;
    auto& gshp = *(glshape*)shp->gl_data;
    clear_glbuffer(gshp.pos);
    clear_glbuffer(gshp.norm);
    clear_glbuffer(gshp.texcoord);
    clear_glbuffer(gshp.color);
    clear_glbuffer(gshp.tangsp);
    clear_glbuffer(gshp.points);
    clear_glbuffer(gshp.lines);
    clear_glbuffer(gshp.triangles);
    clear_glbuffer(gshp.quads);
    clear_glbuffer(gshp.beziers);
    clear_glbuffer(gshp.edges);
    delete (glshape*)shp->gl_data;
}
void clear_gldata(const std::shared_ptr<scene>& scn) {
    for (auto txt : scn->textures) clear_gldata(txt);
    for (auto shp : scn->shapes) clear_gldata(shp);
}

// Draw a shape
void draw_glshape(const std::shared_ptr<shape>& shp,
    const std::shared_ptr<material>& mat, const mat4f& xform, bool highlighted,
    const glsurface_program& prog, bool eyelight, bool wireframe, bool edges) {
    enable_glwireframe(wireframe);
    begin_glsurface_shape(prog, xform);

    auto txt = [](const std::shared_ptr<texture>& txt) -> gltexture_info {
        auto ginfo = gltexture_info();
        if (!txt) return ginfo;
        ginfo.txt = *(gltexture*)txt->gl_data;
        return ginfo;
    };

    auto& gshp = *(glshape*)shp->gl_data;
    auto faceted = shp->norm.empty();

    set_glsurface_material(prog, mat->ke, mat->kd, mat->ks, mat->rs, mat->op,
        txt(mat->ke_txt), txt(mat->kd_txt), txt(mat->ks_txt), txt(mat->rs_txt),
        txt(mat->op_txt), txt(mat->norm_txt), mat->double_sided,
        mat->base_metallic, mat->gltf_textures);

    set_glsurface_vert(
        prog, gshp.pos, gshp.norm, gshp.texcoord, gshp.color, gshp.tangsp);

    if (!is_glbuffer_empty(gshp.points)) {
        set_glsurface_elems(prog, glelem_type::point, faceted);
        draw_glelems(gshp.points);
    }
    if (!is_glbuffer_empty(gshp.lines)) {
        set_glsurface_elems(prog, glelem_type::line, faceted);
        draw_glelems(gshp.lines);
    }
    if (!is_glbuffer_empty(gshp.triangles)) {
        set_glsurface_elems(prog, glelem_type::triangle, faceted);
        draw_glelems(gshp.triangles);
    }
    if (!is_glbuffer_empty(gshp.quads)) {
        set_glsurface_elems(prog, glelem_type::triangle, faceted);
        draw_glelems(gshp.quads);
    }
    if (!is_glbuffer_empty(gshp.beziers)) {
        set_glsurface_elems(prog, glelem_type::line, faceted);
        draw_glelems(gshp.beziers);
    }

    if ((edges && !wireframe) || highlighted) {
        enable_glculling(false);
        set_glsurface_constmaterial(prog,
            (highlighted) ? vec3f{1, 1, 0} : vec3f{0, 0, 0},
            (highlighted) ? 1 : mat->op);
        set_glsurface_normaloffset(prog, 0.01f);
        if (!is_glbuffer_empty(gshp.edges)) draw_glelems(gshp.edges);
    }

    end_glsurface_shape(prog);
    enable_glwireframe(false);
}

// Display a scene
void draw_glscene(const std::shared_ptr<scene>& scn,
    const std::shared_ptr<camera>& cam, const glsurface_program& prog,
    const vec2i& viewport_size, const std::shared_ptr<void>& highlighted,
    bool eyelight, bool wireframe, bool edges, float exposure, float gamma) {
    // begin frame
    enable_gldepth_test(true);
    set_glviewport(viewport_size);

    auto camera_view = frame_to_mat(inverse(cam->frame));
    auto camera_proj =
        (cam->far >= flt_max) ?
            perspective_mat(eval_camera_fovy(cam),
                (float)viewport_size.x / (float)viewport_size.y, cam->near) :
            perspective_mat(eval_camera_fovy(cam),
                (float)viewport_size.x / (float)viewport_size.y, cam->near,
                cam->far);

    begin_glsurface_frame(prog, cam->frame.o, camera_view, camera_proj,
        eyelight, exposure, gamma);

    if (!eyelight) {
        auto lights_pos = std::vector<ygl::vec3f>();
        auto lights_ke = std::vector<ygl::vec3f>();
        auto lights_type = std::vector<int>();
        for (auto lgt : scn->lights) {
            if (lights_pos.size() >= 16) break;
            auto shp = lgt->shp;
            auto bbox = shp->bbox;
            auto pos = (bbox.max + bbox.min) / 2;
            auto area = 0.0f;
            if (!shp->triangles.empty()) {
                for (auto t : shp->triangles)
                    area += triangle_area(
                        shp->pos[t.x], shp->pos[t.y], shp->pos[t.z]);
            } else if (!shp->lines.empty()) {
                for (auto l : shp->lines)
                    area += line_length(shp->pos[l.x], shp->pos[l.y]);
            } else {
                area += shp->pos.size();
            }
            auto ke = lgt->mat->ke * area;
            lights_pos.push_back(transform_point(lgt->frame, pos));
            lights_ke.push_back(ke);
            lights_type.push_back(0);
        }
        set_glsurface_lights(
            prog, {0, 0, 0}, lights_pos, lights_ke, lights_type);
    }

    for (auto ist : scn->instances) {
        draw_glshape(ist->shp, ist->mat, frame_to_mat(ist->frame),
            ist == highlighted || ist->shp == highlighted ||
                ist->mat == highlighted,
            prog, eyelight, wireframe, edges);
    }

    end_glsurface_frame(prog);
    enable_glwireframe(false);
}

float handle_glcamera_turntable_dist(const std::shared_ptr<glwindow>& win,
    const std::shared_ptr<camera>& cam, const std::shared_ptr<scene>& scn,
    const scene_selection& sel) {
    static auto mouse_button_last = false;
    static auto mouse_dist_last = 0.0f;
    auto mouse_button = get_glwindow_mouse_button(win);
    auto alt_down = get_glwindow_alt_key(win);
    // auto shift_down = get_glwindow_shift_key(win);
    // auto ctrl_down = get_glwindow_ctrl_key(win);
    auto widgets_active = get_glwidgets_active(win);

    if (widgets_active || !mouse_button || !alt_down) {
        mouse_button_last = false;
        return 0;
    }

    if (!mouse_button_last) {
        auto center = zero3f;
        if (sel.as<instance>()) {
            auto ist = sel.as<instance>();
            center = transform_point(
                ist->frame, (ist->shp->bbox.min + ist->shp->bbox.max) / 2);
        } else {
            center = (scn->bbox.min + scn->bbox.max) / 2;
        }
        mouse_dist_last = max(0.01f, -dot(center - cam->frame.o, cam->frame.z));
    }
    mouse_button_last = mouse_button;
    return mouse_dist_last;
}

// Handles camera navigation
bool handle_glcamera_turntable(const std::shared_ptr<glwindow>& win,
    const std::shared_ptr<camera>& cam, float dist) {
    static auto mouse_last = zero2f;
    auto mouse_pos = get_glwidnow_mouse_posf(win);
    auto mouse_button = get_glwindow_mouse_button(win);
    auto alt_down = get_glwindow_alt_key(win);
    auto shift_down = get_glwindow_shift_key(win);
    // auto ctrl_down = get_glwindow_ctrl_key(win);
    auto widgets_active = get_glwidgets_active(win);

    if (widgets_active || !alt_down) return false;

    // updated
    auto updated = false;

    // handle mouse and keyboard for navigation
    if (mouse_button) {
        if (mouse_button == 1 && shift_down) mouse_button = 3;
        auto dolly = 0.0f;
        auto pan = zero2f;
        auto rotate = zero2f;
        switch (mouse_button) {
            case 1: rotate = (mouse_pos - mouse_last) / 100.0f; break;
            case 2: dolly = (mouse_pos.x - mouse_last.x) / 100.0f; break;
            case 3: pan = (mouse_pos - mouse_last) / 100.0f; break;
            default: break;
        }
        camera_turntable(cam->frame, dist, rotate, dolly, pan);
        updated = true;
    }

    // record mouse position
    mouse_last = mouse_pos;

    // done
    return updated;
}

// Handles camera navigation
bool handle_glcamera_fps(
    const std::shared_ptr<glwindow>& win, const std::shared_ptr<camera>& cam) {
    static auto mouse_last = zero2f;
    auto mouse_pos = get_glwidnow_mouse_posf(win);
    auto mouse_button = get_glwindow_mouse_button(win);
    auto alt_down = get_glwindow_alt_key(win);
    auto shift_down = get_glwindow_shift_key(win);
    // auto ctrl_down = get_glwindow_ctrl_key(win);
    auto widgets_active = get_glwidgets_active(win);

    if (widgets_active || !alt_down) return false;

    // updated
    auto updated = false;

    // handle mouse and keyboard for navigation
    if (mouse_button) {
        if (mouse_button == 1 && shift_down) mouse_button = 3;
        auto dolly = 0.0f;
        auto pan = zero2f;
        auto rotate = zero2f;
        switch (mouse_button) {
            case 1: rotate = (mouse_pos - mouse_last) / 100.0f; break;
            case 2: dolly = (mouse_pos.x - mouse_last.x) / 100.0f; break;
            case 3: pan = (mouse_pos - mouse_last) / 100.0f; break;
            default: break;
        }
        camera_fps(cam->frame, {0, 0, 0}, rotate);
        updated = true;
    }

    // handle keytboard for navigation
    auto transl = zero3f;
    if (get_glwindow_glkey(win, 'a')) transl.x -= 1;
    if (get_glwindow_glkey(win, 'd')) transl.x += 1;
    if (get_glwindow_glkey(win, 's')) transl.z += 1;
    if (get_glwindow_glkey(win, 'w')) transl.z -= 1;
    if (get_glwindow_glkey(win, 'e')) transl.y += 1;
    if (get_glwindow_glkey(win, 'q')) transl.y -= 1;
    if (transl != zero3f) {
        camera_fps(cam->frame, transl, {0, 0});
        updated = true;
    }

    // record mouse position
    mouse_last = mouse_pos;

    // done
    return updated;
}

// Handle scene selection
bool handle_glscene_selection(const std::shared_ptr<glwindow>& win,
    const std::shared_ptr<scene>& scn, const std::shared_ptr<camera>& cam,
    const vec2i& imsize, const frame2f& imframe, scene_selection& sel) {
    auto mouse_pos = get_glwidnow_mouse_posf(win);
    auto mouse_button = get_glwindow_mouse_button(win);
    auto mouse_mods = get_glwindow_alt_key(win) || get_glwindow_ctrl_key(win) ||
                      get_glwindow_shift_key(win);
    auto widgets_active = get_glwidgets_active(win);

    if (mouse_button != 1 || widgets_active || mouse_mods) return false;
    auto ij = get_glimage_coords(mouse_pos, imframe, imsize);
    if (ij.x < 0 || ij.x >= imsize.x || ij.y < 0 || ij.y >= imsize.y)
        return false;
    auto ray = eval_camera_ray(cam, ij, imsize, {0.5f, 0.5f}, zero2f);
    auto isec = intersect_ray(scn, ray);
    if (!isec.ist) return false;
    sel = isec.ist;
    return true;
}

static const std::unordered_map<std::string, vec4f>
    draw_visitor_highlight_colors = {{"red", {1, 0.5f, 0.5f, 1}},
        {"green", {0.5f, 1, 0.5f, 1}}, {"blue", {0.5f, 0.5f, 1, 1}}};

vec4f get_highlight_color(
    const std::unordered_map<std::string, std::string>& highlights,
    const std::string& name) {
    if (!highlights.empty() && highlights.find(name) != highlights.end()) {
        return draw_visitor_highlight_colors.at(highlights.at(name));
    }
    return zero4f;
}

template <typename T>
void draw_scene_tree_glwidgets_rec(const std::string& lbl_, T* val,
    scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {}

template <typename T>
void draw_glwidgets_scene_tree(const std::string& lbl_, T* val,
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
void draw_scene_tree_glwidgets_rec(const std::string& lbl_,
    const std::shared_ptr<T>& val, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {}

template <typename T>
void draw_glwidgets_scene_tree(const std::string& lbl_,
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
void draw_scene_tree_glwidgets_rec<instance>(const std::string& lbl_,
    const std::shared_ptr<instance>& val, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    draw_glwidgets_scene_tree("shp", val->shp, sel, highlights);
    draw_glwidgets_scene_tree("sbd", val->sbd, sel, highlights);
    draw_glwidgets_scene_tree("mat", val->mat, sel, highlights);
}

template <>
void draw_scene_tree_glwidgets_rec<material>(const std::string& lbl_,
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
void draw_scene_tree_glwidgets_rec<environment>(const std::string& lbl_,
    const std::shared_ptr<environment>& val, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    draw_glwidgets_scene_tree("ke", val->ke_txt, sel, highlights);
}
template <>
void draw_scene_tree_glwidgets_rec<node>(const std::string& lbl_,
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
void draw_scene_tree_glwidgets_rec<animation>(const std::string& lbl_,
    const std::shared_ptr<animation>& val, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    auto tid = 0;
    for (auto tg : val->targets) {
        draw_glwidgets_scene_tree(
            "tg" + std::to_string(tid++), tg, sel, highlights);
    }
}

void draw_glwidgets_scene_tree(const std::shared_ptr<scene>& scn,
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
bool draw_glwidgets_scene_inspector(
    const std::shared_ptr<camera>& val, const std::shared_ptr<scene>& scn) {
    auto edited = 0;
    edited += ImGui::InputText("name", &val->name);
    edited += ImGui::DragFloat3("frame.x", &val->frame.x.x, -1, 1);
    edited += ImGui::DragFloat3("frame.y", &val->frame.y.x, -1, 1);
    edited += ImGui::DragFloat3("frame.z", &val->frame.z.x, -1, 1);
    edited += ImGui::DragFloat3("frame.o", &val->frame.o.x, -10, 10);
    edited += ImGui::Checkbox("ortho", &val->ortho);
    edited += ImGui::DragFloat("width", &val->width, 0.01f, 1);
    edited += ImGui::DragFloat("height", &val->height, 0.01f, 1);
    edited += ImGui::DragFloat("focal", &val->focal, 0.01f, 1);
    edited += ImGui::DragFloat("focus", &val->focus, 0.01f, 1000);
    edited += ImGui::DragFloat("aperture", &val->aperture, 0, 5);
    edited += ImGui::DragFloat("near", &val->near, 0.01f, 10);
    edited += ImGui::DragFloat("far", &val->far, 10, 10000);
    return edited;
}

/// Visit struct elements.
bool draw_glwidgets_scene_inspector(
    const std::shared_ptr<texture>& val, const std::shared_ptr<scene>& scn) {
    auto edited = 0;
    edited += ImGui::InputText("name", &val->name);
    edited += ImGui::InputText("path", &val->path);
    edited += ImGui::Checkbox("clamp", &val->clamp);
    edited += ImGui::DragFloat("scale", &val->scale);
    return edited;
}

bool draw_glwidgets_scene_inspector(
    const std::shared_ptr<material>& val, const std::shared_ptr<scene>& scn) {
    auto edited = 0;
    edited += ImGui::InputText("name", &val->name);
    edited += ImGui::ColorEdit3("ke", &val->ke.x);  // TODO: HDR
    edited += ImGui::ColorEdit3("kd", &val->kd.x);
    edited += ImGui::ColorEdit3("ks", &val->ks.x);
    edited += ImGui::ColorEdit3("kt", &val->kt.x);
    edited += ImGui::DragFloat("rs", &val->rs);
    edited += ImGui::DragFloat("op", &val->op);
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

bool draw_glwidgets_scene_inspector(
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

bool draw_glwidgets_scene_inspector(
    const std::shared_ptr<subdiv>& val, const std::shared_ptr<scene>& scn) {
    auto edited = 0;
    edited += ImGui::InputText("name", &val->name);
    edited += ImGui::DragInt("level", &val->level, 0, 10);
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

bool draw_glwidgets_scene_inspector(
    const std::shared_ptr<instance>& val, const std::shared_ptr<scene>& scn) {
    auto edited = 0;
    edited += ImGui::InputText("name", &val->name);
    edited += ImGui::DragFloat3("frame.x", &val->frame.x.x, -1, 1);
    edited += ImGui::DragFloat3("frame.y", &val->frame.y.x, -1, 1);
    edited += ImGui::DragFloat3("frame.z", &val->frame.z.x, -1, 1);
    edited += ImGui::DragFloat3("frame.o", &val->frame.o.x, -10, 10);
    edited += ImGui::Combo("shp", &val->shp, scn->shapes, true);
    edited += ImGui::Combo("sbd", &val->sbd, scn->subdivs, true);
    edited += ImGui::Combo("mat", &val->mat, scn->materials, true);
    return edited;
}

bool draw_glwidgets_scene_inspector(const std::shared_ptr<environment>& val,
    const std::shared_ptr<scene>& scn) {
    auto edited = 0;
    edited += ImGui::InputText("name", &val->name);
    edited += ImGui::DragFloat3("frame.x", &val->frame.x.x, -1, 1);
    edited += ImGui::DragFloat3("frame.y", &val->frame.y.x, -1, 1);
    edited += ImGui::DragFloat3("frame.z", &val->frame.z.x, -1, 1);
    edited += ImGui::DragFloat3("frame.o", &val->frame.o.x, -10, 10);
    edited += ImGui::ColorEdit4("ke", &val->ke.x);  // TODO: HDR
    edited += ImGui::Combo("ke txt", &val->ke_txt, scn->textures, true);
    return edited;
}

bool draw_glwidgets_scene_inspector(
    const std::shared_ptr<node>& val, const std::shared_ptr<scene>& scn) {
    auto edited = 0;
    edited += ImGui::InputText("name", &val->name);
    edited += ImGui::Combo("parent", &val->parent, scn->nodes, true);
    edited += ImGui::DragFloat3("frame.x", &val->frame.x.x, -1, 1);
    edited += ImGui::DragFloat3("frame.y", &val->frame.y.x, -1, 1);
    edited += ImGui::DragFloat3("frame.z", &val->frame.z.x, -1, 1);
    edited += ImGui::DragFloat3("frame.o", &val->frame.o.x, -10, 10);
    edited += ImGui::DragFloat3("translation", &val->translation.x);
    edited += ImGui::DragFloat4("rotation", &val->rotation.x, -1, 1);
    edited += ImGui::DragFloat3("scale", &val->scale.x, 0, 10);
    edited += ImGui::Combo("cam", &val->cam, scn->cameras, true);
    edited += ImGui::Combo("ist", &val->ist, scn->instances, true);
    edited += ImGui::Combo("env", &val->env, scn->environments, true);
    return edited;
}

bool draw_glwidgets_scene_inspector(
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

bool draw_glwidgets_scene_tree(const std::string& lbl,
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

bool draw_glwidgets_scene_inspector(const std::string& lbl,
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
