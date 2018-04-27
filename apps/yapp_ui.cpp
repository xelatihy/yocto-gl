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
#include "../yocto/yocto_bvh.h"
#include "../yocto/yocto_shape.h"

namespace ygl {

static const std::map<material_type, std::string>& material_type_names() {
    static auto names = std::map<material_type, std::string>{
        {material_type::specular_roughness, "specular_roughness"},
        {material_type::metallic_roughness, "metallic_roughness"},
        {material_type::specular_glossiness, "specular_glossiness"},
    };
    return names;
}

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

void update_gldata(texture* txt) {
    if (!txt->gl_data) txt->gl_data = new gltexture();
    auto& gtxt = *(gltexture*)txt->gl_data;
    if (!txt->hdr.empty()) {
        update_gltexture(
            gtxt, txt->width, txt->height, txt->hdr, true, true, true);
    }
    if (!txt->ldr.empty()) {
        update_gltexture(
            gtxt, txt->width, txt->height, txt->ldr, true, true, true);
    }
}
void update_gldata(shape* shp) {
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
    if (!shp->quads_pos.empty()) {
        auto pos = std::vector<vec3f>();
        auto norm = std::vector<vec3f>();
        auto texcoord = std::vector<vec2f>();
        auto quads = std::vector<vec4i>();
        std::tie(quads, pos, norm, texcoord) =
            convert_face_varying(shp->quads_pos, shp->quads_norm,
                shp->quads_texcoord, shp->pos, shp->norm, shp->texcoord);
        update_vert_buffer(gshp.pos, pos);
        update_vert_buffer(gshp.norm, norm);
        update_vert_buffer(gshp.texcoord, texcoord);
        update_vert_buffer(gshp.color, std::vector<vec4f>{});
        update_vert_buffer(gshp.tangsp, std::vector<vec4f>{});
        update_elem_buffer(gshp.quads, convert_quads_to_triangles(quads));
        update_elem_buffer(gshp.edges, get_edges({}, {}, quads));
    } else {
        update_vert_buffer(gshp.pos, shp->pos);
        update_vert_buffer(gshp.norm, shp->norm);
        update_vert_buffer(gshp.texcoord, shp->texcoord);
        update_vert_buffer(gshp.color, shp->color);
        update_vert_buffer(gshp.tangsp, shp->tangsp);
        update_elem_buffer(gshp.points, shp->points);
        update_elem_buffer(gshp.lines, shp->lines);
        update_elem_buffer(gshp.triangles, shp->triangles);
        update_elem_buffer(gshp.quads, convert_quads_to_triangles(shp->quads));
        update_elem_buffer(gshp.beziers, convert_bezier_to_lines(shp->beziers));
        update_elem_buffer(
            gshp.edges, get_edges({}, shp->triangles, shp->quads));
    }
}
void update_gldata(scene* scn) {
    for (auto txt : scn->textures) update_gldata(txt);
    for (auto shp : scn->shapes) update_gldata(shp);
}
void clear_gldata(texture* txt) {
    if (!txt->gl_data) return;
    clear_gltexture(*(gltexture*)txt->gl_data);
    delete (gltexture*)txt->gl_data;
}
void clear_gldata(shape* shp) {
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
void clear_gldata(scene* scn) {
    for (auto txt : scn->textures) clear_gldata(txt);
    for (auto shp : scn->shapes) clear_gldata(shp);
}

// Draw a shape
void draw_glshape(const shape* shp, const material* mat, const mat4f& xform,
    bool highlighted, const glsurface_program& prog, bool eyelight,
    bool wireframe, bool edges, bool cutout, bool cull) {
    static auto default_material = material();
    default_material.kd = {0.2f, 0.2f, 0.2f};
    static auto mtypes = std::unordered_map<material_type, int>{
        {material_type::specular_roughness, 1},
        {material_type::metallic_roughness, 2},
        {material_type::specular_glossiness, 3}};

    if (!mat) mat = &default_material;

    enable_glwireframe(wireframe);
    enable_glculling(cull && !mat->double_sided);
    begin_glsurface_shape(prog, xform);

    auto txt = [](const texture_info& info) -> gltexture_info {
        auto ginfo = gltexture_info();
        if (!info.txt) return ginfo;
        ginfo.txt = *(gltexture*)info.txt->gl_data;
        return ginfo;
    };

    auto& gshp = *(glshape*)shp->gl_data;
    auto faceted = shp->norm.empty();

    set_glsurface_material(prog, mtypes.at(mat->type), mat->ke, mat->kd,
        mat->ks, mat->rs, mat->op, txt(mat->ke_txt), txt(mat->kd_txt),
        txt(mat->ks_txt), txt(mat->rs_txt), txt(mat->norm_txt),
        txt(mat->occ_txt), false, mat->double_sided, cutout);

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
        if (is_glbuffer_empty(gshp.edges)) draw_glelems(gshp.edges);
    }

    end_glsurface_shape(prog);
    enable_glculling(false);
    enable_glwireframe(false);
}

// Display a scene
void draw_glscene(const scene* scn, const camera* cam,
    const glsurface_program& prog, const vec2i& viewport_size,
    const void* highlighted, bool eyelight, bool wireframe, bool edges,
    bool cutout, float exposure, float gamma, bool cull) {
    // begin frame
    enable_gldepth_test(true);
    set_glviewport(viewport_size);

    auto camera_xform = frame_to_mat(cam->frame);
    auto camera_view = frame_to_mat(inverse(cam->frame));
    auto camera_proj = perspective_mat(cam->yfov,
        (float)viewport_size.x / (float)viewport_size.y, cam->near, cam->far);

    begin_glsurface_frame(prog, camera_xform, camera_view, camera_proj,
        eyelight, exposure, gamma);

    if (!eyelight) {
        auto lights_pos = std::vector<ygl::vec3f>();
        auto lights_ke = std::vector<ygl::vec3f>();
        auto lights_type = std::vector<int>();
        for (auto lgt : scn->lights) {
            if (lights_pos.size() >= 16) break;
            if (!lgt->ist) continue;
            auto shp = lgt->ist->shp;
            if (!shp->points.empty()) {
                for (auto p : shp->points) {
                    if (lights_pos.size() >= 16) break;
                    lights_pos.push_back(
                        transform_point(lgt->ist->frame, shp->pos[p]));
                    lights_ke.push_back(lgt->ist->mat->ke);
                    lights_type.push_back(0);
                }
            } else {
                auto bbox = shp->bbox;
                auto pos = (bbox.max + bbox.min) / 2;
                auto area = 0.0f;
                for (auto l : shp->lines)
                    area += line_length(shp->pos[l.x], shp->pos[l.y]);
                for (auto t : shp->triangles)
                    area += triangle_area(
                        shp->pos[t.x], shp->pos[t.y], shp->pos[t.z]);
                for (auto t : shp->quads)
                    area += quad_area(shp->pos[t.x], shp->pos[t.y],
                        shp->pos[t.z], shp->pos[t.w]);
                auto ke = lgt->ist->mat->ke * area;
                lights_pos.push_back(transform_point(lgt->ist->frame, pos));
                lights_ke.push_back(ke);
                lights_type.push_back(0);
            }
        }
        set_glsurface_lights(
            prog, {0, 0, 0}, lights_pos, lights_ke, lights_type);
    }

    for (auto ist : scn->instances) {
        draw_glshape(ist->shp, ist->mat, frame_to_mat(ist->frame),
            ist == highlighted || ist->shp == highlighted ||
                ist->mat == highlighted,
            prog, eyelight, wireframe, edges, cutout, cull);
    }

    end_glsurface_frame(prog);
    enable_glwireframe(false);
}

// Handles camera navigation
bool handle_glcamera_navigation(
    glwindow* win, camera* cam, bool navigation_fps) {
    static auto mouse_last = zero2f;
    auto mouse_pos = get_glwidnow_mouse_posf(win);
    auto mouse_button = get_glwindow_mouse_button(win);
    auto alt_down = get_glwindow_alt_key(win);

    // updated
    auto updated = false;

    // handle mouse and keyboard for navigation
    if (mouse_button && alt_down && !get_glwidgets_active(win)) {
        if (mouse_button == 1 && get_glwindow_shift_key(win)) mouse_button = 3;
        if (navigation_fps) {
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
        } else {
            auto dolly = 0.0f;
            auto pan = zero2f;
            auto rotate = zero2f;
            switch (mouse_button) {
                case 1: rotate = (mouse_pos - mouse_last) / 100.0f; break;
                case 2: dolly = (mouse_pos.x - mouse_last.x) / 100.0f; break;
                case 3: pan = (mouse_pos - mouse_last) / 100.0f; break;
                default: break;
            }
            camera_turntable(cam->frame, cam->focus, rotate, dolly, pan);
            updated = true;
        }
    }

    // handle keytboard for navigation
    if (!get_glwidgets_active(win) && navigation_fps) {
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
    }

    // record mouse position
    mouse_last = mouse_pos;

    // done
    return updated;
}

// Handle scene selection
bool handle_glscene_selection(glwindow* win, const scene* scn,
    const camera* cam, int res, const frame2f& imframe, scene_selection& sel) {
    auto mouse_pos = get_glwidnow_mouse_posf(win);
    auto mouse_button = get_glwindow_mouse_button(win);

    if (!(mouse_button == 1 && !get_glwidgets_active(win))) return false;
    auto ij = get_glimage_coords(
        mouse_pos, imframe, {(int)round(cam->aspect * res), res});
    if (ij.x < 0 || ij.x >= (int)round(res * cam->aspect) || ij.y < 0 ||
        ij.y >= res)
        return false;
    auto ray = eval_camera_ray(cam, ij, res, {0.5f, 0.5f}, zero2f);
    auto iid = 0, eid = 0;
    auto ray_t = 0.0f;
    auto euv = zero2f;
    if (!intersect_bvh(scn->bvh, ray, false, ray_t, iid, eid, euv))
        return false;
    if (scn->nodes.empty()) {
        sel = scn->shapes[iid];
    } else {
        sel = scn->nodes[iid];
    }
    return true;
}

// Implementation of camera selection
bool draw_glwidgets_camera_inspector(
    glwindow* win, const std::string& lbl, camera* cam) {
    if (!cam) return false;
    auto edited = 0;
    auto from = cam->frame.o;
    auto to = cam->frame.o - cam->frame.z * cam->focus;
    edited += draw_glwidgets_text(win, lbl + " name", cam->name);
    edited += draw_glwidgets_dragbox(win, lbl + " from", from, -10, 10);
    edited += draw_glwidgets_dragbox(win, lbl + " to", to, -10, 10);
    edited += draw_glwidgets_dragbox(win, lbl + " yfov", cam->yfov, 0.1, 10);
    edited += draw_glwidgets_dragbox(win, lbl + " aspect", cam->aspect, 1, 3);
    edited +=
        draw_glwidgets_dragbox(win, lbl + " aperture", cam->aperture, 0, 1);
    edited += draw_glwidgets_dragbox(win, lbl + " focus", cam->focus, 0.1, 100);
    edited += draw_glwidgets_dragbox(win, lbl + " near", cam->near, 0.01f, 1);
    edited += draw_glwidgets_dragbox(win, lbl + " far", cam->far, 1, 1000);
    if (edited) cam->frame = lookat_frame(from, to, {0, 1, 0});
    return edited;
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
void draw_scene_tree_glwidgets_rec(glwindow* win, const std::string& lbl_,
    T* val, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {}

template <typename T>
void draw_glwidgets_scene_tree(glwindow* win, const std::string& lbl_, T* val,
    scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    if (!val) return;
    auto lbl = val->name;
    if (!lbl_.empty()) lbl = lbl_ + ": " + val->name;
    auto selection = sel.ptr;
    auto color = get_highlight_color(highlights, val->name);
    if (color != zero4f) push_glwidgets_style(win, color);
    auto open = begin_glwidgets_tree(win, lbl, selection, val);
    if (color != zero4f) pop_glwidgets_style(win);
    if (selection == val) sel = val;
    if (open) {
        draw_scene_tree_glwidgets_rec(win, lbl_, val, sel, highlights);
        end_glwidgets_tree(win);
    }
}

template <>
void draw_scene_tree_glwidgets_rec<instance>(glwindow* win,
    const std::string& lbl_, instance* val, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    draw_glwidgets_scene_tree(win, "shp", val->shp, sel, highlights);
    draw_glwidgets_scene_tree(win, "mat", val->mat, sel, highlights);
}

template <>
void draw_scene_tree_glwidgets_rec<material>(glwindow* win,
    const std::string& lbl_, material* val, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    draw_glwidgets_scene_tree(win, "ke", val->ke_txt.txt, sel, highlights);
    draw_glwidgets_scene_tree(win, "kd", val->kd_txt.txt, sel, highlights);
    draw_glwidgets_scene_tree(win, "ks", val->ks_txt.txt, sel, highlights);
    draw_glwidgets_scene_tree(win, "kr", val->kr_txt.txt, sel, highlights);
    draw_glwidgets_scene_tree(win, "rs", val->rs_txt.txt, sel, highlights);
    draw_glwidgets_scene_tree(win, "op", val->op_txt.txt, sel, highlights);
    draw_glwidgets_scene_tree(win, "bump", val->bump_txt.txt, sel, highlights);
    draw_glwidgets_scene_tree(win, "disp", val->disp_txt.txt, sel, highlights);
    draw_glwidgets_scene_tree(win, "norm", val->norm_txt.txt, sel, highlights);
    draw_glwidgets_scene_tree(win, "occ", val->occ_txt.txt, sel, highlights);
}
template <>
void draw_scene_tree_glwidgets_rec<environment>(glwindow* win,
    const std::string& lbl_, environment* val, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    draw_glwidgets_scene_tree(win, "ke", val->ke_txt.txt, sel, highlights);
}
template <>
void draw_scene_tree_glwidgets_rec<node>(glwindow* win, const std::string& lbl_,
    node* val, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    draw_glwidgets_scene_tree(win, "ist", val->ist, sel, highlights);
    draw_glwidgets_scene_tree(win, "cam", val->cam, sel, highlights);
    draw_glwidgets_scene_tree(win, "env", val->env, sel, highlights);
    draw_glwidgets_scene_tree(win, "par", val->parent, sel, highlights);
    auto cid = 0;
    for (auto ch : val->children) {
        draw_glwidgets_scene_tree(
            win, "ch" + std::to_string(cid++), ch, sel, highlights);
    }
}
template <>
void draw_scene_tree_glwidgets_rec<animation>(glwindow* win,
    const std::string& lbl_, animation* val, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    auto tid = 0;
    for (auto tg : val->targets) {
        draw_glwidgets_scene_tree(
            win, "tg" + std::to_string(tid++), tg, sel, highlights);
    }
}

void draw_glwidgets_scene_tree(glwindow* win, const std::string& lbl,
    texture_info& val, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    draw_glwidgets_scene_tree(win, lbl, val.txt, sel, highlights);
}

void draw_glwidgets_scene_tree(glwindow* win, scene* scn, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    if (!scn->cameras.empty() && begin_glwidgets_tree(win, "cameras")) {
        for (auto v : scn->cameras)
            draw_glwidgets_scene_tree(win, "", v, sel, highlights);
        end_glwidgets_tree(win);
    }
    if (!scn->shapes.empty() && begin_glwidgets_tree(win, "shapes")) {
        for (auto v : scn->shapes)
            draw_glwidgets_scene_tree(win, "", v, sel, highlights);
        end_glwidgets_tree(win);
    }
    if (!scn->instances.empty() && begin_glwidgets_tree(win, "instances")) {
        for (auto v : scn->instances)
            draw_glwidgets_scene_tree(win, "", v, sel, highlights);
        end_glwidgets_tree(win);
    }
    if (!scn->materials.empty() && begin_glwidgets_tree(win, "materials")) {
        for (auto v : scn->materials)
            draw_glwidgets_scene_tree(win, "", v, sel, highlights);
        end_glwidgets_tree(win);
    }
    if (!scn->textures.empty() && begin_glwidgets_tree(win, "textures")) {
        for (auto v : scn->textures)
            draw_glwidgets_scene_tree(win, "", v, sel, highlights);
        end_glwidgets_tree(win);
    }
    if (!scn->environments.empty() &&
        begin_glwidgets_tree(win, "environments")) {
        for (auto v : scn->environments)
            draw_glwidgets_scene_tree(win, "", v, sel, highlights);
        end_glwidgets_tree(win);
    }
    if (!scn->nodes.empty() && begin_glwidgets_tree(win, "nodes")) {
        for (auto v : scn->nodes)
            draw_glwidgets_scene_tree(win, "", v, sel, highlights);
        end_glwidgets_tree(win);
    }
    if (!scn->animations.empty() && begin_glwidgets_tree(win, "animations")) {
        for (auto v : scn->animations)
            draw_glwidgets_scene_tree(win, "", v, sel, highlights);
        end_glwidgets_tree(win);
    }
}

template <typename T>
void draw_glwidgets_label(glwindow* win, const std::string& lbl,
    const std::vector<T>& val, bool skip_if_empty = true) {
    if (skip_if_empty && val.empty()) return;
    draw_glwidgets_label(win, lbl, std::to_string(val.size()));
}

/// Visit struct elements.
bool draw_glwidgets_scene_inspector(glwindow* win, camera* val, scene* scn) {
    auto edited = 0;
    edited += draw_glwidgets_text(win, "name", val->name);
    edited += draw_glwidgets_dragbox(win, "frame", val->frame);
    edited += draw_glwidgets_checkbox(win, "ortho", val->ortho);
    edited += draw_glwidgets_dragbox(win, "yfov", val->yfov, 0.1f, 10);
    edited += draw_glwidgets_dragbox(win, "aspect", val->aspect, 1, 3);
    edited += draw_glwidgets_dragbox(win, "focus", val->focus, 0.01f, 1000);
    edited += draw_glwidgets_dragbox(win, "aperture", val->aperture, 0, 5);
    edited += draw_glwidgets_dragbox(win, "near", val->near, 0.01f, 10);
    edited += draw_glwidgets_dragbox(win, "far", val->far, 10, 10000);
    return edited;
}

/// Visit struct elements.
bool draw_glwidgets_scene_inspector(glwindow* win, texture* val, scene* scn) {
    auto edited = 0;
    edited += draw_glwidgets_text(win, "name", val->name);
    edited += draw_glwidgets_text(win, "path", val->path);
    draw_glwidgets_label(win, "ldr", val->ldr);
    draw_glwidgets_label(win, "hdr", val->hdr);
    return edited;
}

bool draw_glwidgets_scene_inspector(
    glwindow* win, const std::string& lbl, texture_info& val, scene* scn) {
    auto edited = 0;
    edited += draw_glwidgets_combobox(
        win, lbl + " txt", val.txt, scn->textures, true);
    edited += draw_glwidgets_checkbox(win, lbl + " wrap_s", val.wrap_s);
    edited += draw_glwidgets_checkbox(win, lbl + " wrap_t", val.wrap_t);
    edited += draw_glwidgets_checkbox(win, lbl + " linear", val.linear);
    edited += draw_glwidgets_checkbox(win, lbl + " mipmap", val.mipmap);
    edited += draw_glwidgets_dragbox(win, lbl + " scale", val.scale);
    return edited;
}

bool draw_glwidgets_scene_inspector(glwindow* win, material* val, scene* scn) {
    auto edited = 0;
    edited += draw_glwidgets_text(win, "name", val->name);
    edited += draw_glwidgets_checkbox(win, "double_sided", val->double_sided);
    edited +=
        draw_glwidgets_combobox(win, "type", val->type, material_type_names());
    edited += draw_hdr_color_widget(win, "ke", val->ke);
    edited += draw_glwidgets_colorbox(win, "kd", val->kd);
    edited += draw_glwidgets_colorbox(win, "ks", val->ks);
    edited += draw_glwidgets_colorbox(win, "kr", val->kr);
    edited += draw_glwidgets_colorbox(win, "kt", val->kt);
    edited += draw_glwidgets_dragbox(win, "rs", val->rs);
    edited += draw_glwidgets_dragbox(win, "op", val->op);
    edited += draw_glwidgets_scene_inspector(win, "ke", val->ke_txt, scn);
    edited += draw_glwidgets_scene_inspector(win, "kd", val->kd_txt, scn);
    edited += draw_glwidgets_scene_inspector(win, "ks", val->ks_txt, scn);
    edited += draw_glwidgets_scene_inspector(win, "kr", val->kr_txt, scn);
    edited += draw_glwidgets_scene_inspector(win, "kt", val->kt_txt, scn);
    edited += draw_glwidgets_scene_inspector(win, "rs", val->rs_txt, scn);
    edited += draw_glwidgets_scene_inspector(win, "op", val->op_txt, scn);
    edited += draw_glwidgets_scene_inspector(win, "bump", val->bump_txt, scn);
    edited += draw_glwidgets_scene_inspector(win, "disp", val->disp_txt, scn);
    edited += draw_glwidgets_scene_inspector(win, "norm", val->norm_txt, scn);
    edited += draw_glwidgets_scene_inspector(win, "occ", val->occ_txt, scn);
    return edited;
}

bool draw_glwidgets_scene_inspector(glwindow* win, shape* val, scene* scn) {
    auto edited = 0;
    edited += draw_glwidgets_text(win, "name", val->name);
    edited += draw_glwidgets_text(win, "path", val->path);
    draw_glwidgets_label(win, "points", val->points);
    draw_glwidgets_label(win, "lines", val->lines);
    draw_glwidgets_label(win, "triangles", val->triangles);
    draw_glwidgets_label(win, "quads", val->quads);
    draw_glwidgets_label(win, "quads_pos", val->quads_pos);
    draw_glwidgets_label(win, "quads_norm", val->quads_norm);
    draw_glwidgets_label(win, "quads_texcoord", val->quads_texcoord);
    draw_glwidgets_label(win, "beziers", val->beziers);
    draw_glwidgets_label(win, "pos", val->pos);
    draw_glwidgets_label(win, "norm", val->norm);
    draw_glwidgets_label(win, "texcoord", val->texcoord);
    draw_glwidgets_label(win, "texcoord1", val->texcoord1);
    draw_glwidgets_label(win, "color", val->color);
    draw_glwidgets_label(win, "radius", val->radius);
    draw_glwidgets_label(win, "tangsp", val->tangsp);
    draw_glwidgets_label(win, "subdivision", std::to_string(val->subdivision));
    draw_glwidgets_label(
        win, "catmullclark", (val->catmullclark) ? "true" : "false");
    return edited;
}

bool draw_glwidgets_scene_inspector(glwindow* win, instance* val, scene* scn) {
    auto edited = 0;
    edited += draw_glwidgets_text(win, "name", val->name);
    edited += draw_glwidgets_dragbox(win, "frame", val->frame);
    edited += draw_glwidgets_combobox(win, "shp", val->shp, scn->shapes, true);
    edited +=
        draw_glwidgets_combobox(win, "mat", val->mat, scn->materials, true);
    return edited;
}

bool draw_glwidgets_scene_inspector(
    glwindow* win, environment* val, scene* scn) {
    auto edited = 0;
    edited += draw_glwidgets_text(win, "name", val->name);
    edited += draw_glwidgets_dragbox(win, "frame", val->frame);
    edited += draw_hdr_color_widget(win, "ke", val->ke);
    edited += draw_glwidgets_scene_inspector(win, "ke", val->ke_txt, scn);
    return edited;
}

bool draw_glwidgets_scene_inspector(glwindow* win, node* val, scene* scn) {
    auto edited = 0;
    edited += draw_glwidgets_text(win, "name", val->name);
    edited +=
        draw_glwidgets_combobox(win, "parent", val->parent, scn->nodes, true);
    edited += draw_glwidgets_dragbox(win, "frame", val->frame);
    edited += draw_glwidgets_dragbox(win, "translation", val->translation);
    edited += draw_glwidgets_dragbox(win, "rotation", val->rotation, -1, 1);
    edited += draw_glwidgets_dragbox(win, "scale", val->scale, 0, 10);
    edited += draw_glwidgets_combobox(win, "cam", val->cam, scn->cameras, true);
    edited +=
        draw_glwidgets_combobox(win, "ist", val->ist, scn->instances, true);
    edited +=
        draw_glwidgets_combobox(win, "env", val->env, scn->environments, true);
    return edited;
}

bool draw_glwidgets_scene_inspector(glwindow* win, animation* val, scene* scn) {
    auto edited = 0;
    edited += draw_glwidgets_text(win, "name", val->name);
    edited += draw_glwidgets_text(win, "path", val->path);
    edited += draw_glwidgets_text(win, "group", val->group);
    edited +=
        draw_glwidgets_combobox(win, "type", val->type, animation_type_names());
    draw_glwidgets_label(win, "times", val->times);
    draw_glwidgets_label(win, "translation", val->translation);
    draw_glwidgets_label(win, "rotation", val->rotation);
    draw_glwidgets_label(win, "scale", val->scale);
    draw_glwidgets_label(win, "weights", val->weights);
    draw_glwidgets_label(win, "targets", val->targets);
    return edited;
}

bool draw_glwidgets_scene_tree(glwindow* win, const std::string& lbl,
    scene* scn, scene_selection& sel,
    std::vector<ygl::scene_selection>& update_list,
    const std::unordered_map<std::string, std::string>& inspector_highlights) {
    if (!scn) return false;
    push_glwidgets_groupid(win, scn);
    // begin_glwidgets_scrollarea(win, "scene #$%^!@", 240, false);
    draw_glwidgets_scene_tree(win, scn, sel, inspector_highlights);
    // end_glwidgets_scrollarea(win);

    auto update_len = update_list.size();
#if 0
    if (test_scn) {
        draw_add_elem_glwidgets(
            win, scn, "cam", scn->cameras, test_scn->cameras, sel, update_list);
        draw_add_elem_glwidgets(win, scn, "txt", scn->textures,
            test_scn->textures, sel, update_list);
        draw_add_elem_glwidgets(win, scn, "mat", scn->materials,
            test_scn->materials, sel, update_list);
        draw_add_elem_glwidgets(
            win, scn, "shp", scn->shapes, test_scn->shapes, sel, update_list);
        draw_add_elem_glwidgets(win, scn, "ist", scn->instances,
            test_scn->instances, sel, update_list);
        draw_add_elem_glwidgets(
            win, scn, "nde", scn->nodes, test_scn->nodes, sel, update_list);
        draw_add_elem_glwidgets(win, scn, "env", scn->environments,
            test_scn->environments, sel, update_list);
        draw_add_elem_glwidgets(win, scn, "anim", scn->animations,
            test_scn->animations, sel, update_list);
    }
#endif

    pop_glwidgets_groupid(win);
    return update_list.size() != update_len;
}

bool draw_glwidgets_scene_inspector(glwindow* win, const std::string& lbl,
    scene* scn, scene_selection& sel,
    std::vector<ygl::scene_selection>& update_list,
    const std::unordered_map<std::string, std::string>& inspector_highlights) {
    if (!scn || !sel.ptr) return false;
    push_glwidgets_groupid(win, sel.ptr);

    auto update_len = update_list.size();

    auto edited = false;
    if (sel.as<camera>())
        edited = draw_glwidgets_scene_inspector(win, sel.as<camera>(), scn);
    if (sel.as<shape>())
        edited = draw_glwidgets_scene_inspector(win, sel.as<shape>(), scn);
    if (sel.as<texture>())
        edited = draw_glwidgets_scene_inspector(win, sel.as<texture>(), scn);
    if (sel.as<material>())
        edited = draw_glwidgets_scene_inspector(win, sel.as<material>(), scn);
    if (sel.as<environment>())
        edited =
            draw_glwidgets_scene_inspector(win, sel.as<environment>(), scn);
    if (sel.as<instance>())
        edited = draw_glwidgets_scene_inspector(win, sel.as<instance>(), scn);
    if (sel.as<node>())
        edited = draw_glwidgets_scene_inspector(win, sel.as<node>(), scn);
    if (sel.as<animation>())
        edited = draw_glwidgets_scene_inspector(win, sel.as<animation>(), scn);
    if (edited) update_list.push_back(sel);

    pop_glwidgets_groupid(win);
    return update_list.size() != update_len;
}

}  // namespace ygl
