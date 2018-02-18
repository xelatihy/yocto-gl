//
// LICENSE:
//
// Copyright (c) 2016 -- 2017 Fabio Pellacini
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

#include "../yocto/yocto_gl.h"
#include "../yocto/yocto_gltf.h"
using ygl::log_info;

// OpenGL shape vbo
struct shape_vbo {
    ygl::gl_vertex_buffer pos = {};
    ygl::gl_vertex_buffer norm = {};
    ygl::gl_vertex_buffer texcoord = {};
    ygl::gl_vertex_buffer texcoord1 = {};
    ygl::gl_vertex_buffer color = {};
    ygl::gl_vertex_buffer tangsp = {};
    ygl::gl_vertex_buffer skin_joints = {};
    ygl::gl_vertex_buffer skin_weights = {};
    ygl::gl_element_buffer points = {};
    ygl::gl_element_buffer lines = {};
    ygl::gl_element_buffer triangles = {};
};

// OpenGL state
struct shade_state {
    ygl::gl_stdsurface_program prog = {};
    std::map<void*, ygl::gl_texture> txt;
    std::map<void*, shape_vbo> vbo;
    ygl::gl_lights lights;
};

// Application state
struct app_state {
    // scene data
    ygl::gltf_scene_group* gscn = nullptr;

    // view camera
    ygl::camera* view_cam = nullptr;

    // view/scene selection
    ygl::gltf_node* gcam = nullptr;

    // time and animation
    float time = 0;
    ygl::vec2f time_range = {0, 0};
    bool animate = false;

    // filenames
    std::string filename;
    std::string imfilename;
    std::string outfilename;

    // load options
    bool split_shapes = false;
    bool add_spec_gloss = false;

    // render
    int resolution = 0;
    float exposure = 0, gamma = 2.2f;
    bool filmic = false;
    ygl::vec4f background = {0, 0, 0, 0};

    // lighting
    bool camera_lights = false;
    ygl::vec3f amb = {0, 0, 0};

    // ui
    bool interactive = false;
    bool scene_updated = false;

    // shade
    bool wireframe = false, edges = false;
    bool alpha_cutout = true;
    shade_state* shstate = nullptr;

    // navigation
    bool navigation_fps = false;

    // editing support
    void* selection = nullptr;

    ~app_state() {
        if (shstate) delete shstate;
        if (view_cam) delete view_cam;
    }
};

// loads a scene either OBJ or GLTF
inline bool load_scene(
    app_state* scn, const std::string& filename, bool instances, bool radius) {
    // load scene
    auto ext = ygl::path_extension(filename);
    try {
        scn->gscn = ygl::load_scenes(filename, true, true);
    } catch (const std::exception& e) {
        ygl::log_error(
            "cannot load scene {} with error {}", filename, e.what());
    }
    ygl::add_normals(scn->gscn);
    ygl::add_texture_data(scn->gscn);
    if (scn->add_spec_gloss) ygl::add_spec_gloss(scn->gscn);
    if (radius) ygl::add_radius(scn->gscn, 0.001f);
    ygl::add_tangent_space(scn->gscn);
    if (instances) {
        ygl::add_nodes(scn->gscn);
        ygl::add_scene(scn->gscn);
    }
    ygl::add_names(scn->gscn);
    scn->time_range = ygl::get_animation_bounds(scn->gscn);
    if (!scn->gscn->default_scene && !scn->gscn->scenes.empty()) {
        scn->gscn->default_scene = scn->gscn->scenes[0];
    }
    if (!ygl::get_camera_nodes(scn->gscn->default_scene).empty()) {
        auto cam = ygl::get_camera_nodes(scn->gscn->default_scene)[0];
        scn->view_cam = new ygl::camera();
        scn->view_cam->frame = ygl::mat_to_frame(ygl::mat4f(cam->xform()));
        scn->view_cam->yfov = cam->cam->yfov;
        scn->view_cam->aspect = cam->cam->aspect;
        scn->view_cam->near = cam->cam->near;
        scn->view_cam->far = cam->cam->far;
        scn->view_cam->focus = cam->cam->focus;
        scn->view_cam->aperture = cam->cam->aperture;
    }

#if 1
    if (!scn->view_cam) {
        auto bbox = ygl::bbox3f{ygl::compute_scene_bounds(scn->gscn)};
        auto center = (bbox.max + bbox.min) / 2.0f;
        auto bbox_size = bbox.max - bbox.min;
        auto bbox_msize = length(bbox_size) / 2;
        // auto bbox_msize =
        //     ygl::max(bbox_size[0], ygl::max(bbox_size[1], bbox_size[2]));
        scn->view_cam = new ygl::camera();
        auto camera_dir = ygl::vec3f{1.5f, 0.8f, 1.5f};
        auto from = camera_dir * bbox_msize + center;
        auto to = center;
        auto up = ygl::vec3f{0, 1, 0};
        scn->view_cam->frame = ygl::lookat_frame(from, to, up);
        scn->view_cam->aspect = 16.0f / 9.0f;
        scn->view_cam->yfov = 2 * atanf(0.5f);
        scn->view_cam->aperture = 0;
        scn->view_cam->focus = ygl::length(to - from);
    }
#else
    if (!scn->view_cam) {
        auto bbox = (scn->oscn) ?
                        ygl::bbox3f{yobj::compute_scene_bounds(scn->oscn)} :
                        ygl::bbox3f{ygl::compute_scene_bounds(scn->gscn)};
        auto center = ygl::center(bbox);
        auto bbox_size = ygl::diagonal(bbox);
        scn->view_cam = new ycamera();
        auto from =
            center + ygl::vec3f{0, 0, 1} *
                         (bbox_size.z + ygl::max(bbox_size.x, bbox_size.y));
        auto to = center;
        auto up = ygl::vec3f{0, 1, 0};
        scn->view_cam->frame = ygl::lookat_frame3(from, to, up);
        scn->view_cam->aspect = 16.0f / 9.0f;
        scn->view_cam->yfov = 2 * atanf(0.5f);
        scn->view_cam->aperture = 0;
        scn->view_cam->focus = ygl::length(to - from);
    }
#endif
    return true;
}

// Init shading
inline void update_shade_lights(
    shade_state* st, const ygl::gltf_scene_group* scn) {
    st->lights.pos.clear();
    st->lights.ke.clear();
    st->lights.type.clear();

    auto instances = ygl::get_mesh_nodes(scn->default_scene);

    if (!instances.empty()) {
        for (auto ist : instances) {
            for (auto shp : ist->msh->shapes) {
                if (!shp->mat) continue;
                auto mat = shp->mat;
                if (mat->emission == ygl::zero3f) continue;
                for (auto p : shp->points) {
                    if (st->lights.pos.size() >= 16) break;
                    st->lights.pos.push_back(ygl::transform_point(
                        ygl::mat4f(ist->xform()), shp->pos[p]));
                    st->lights.ke.push_back(mat->emission);
                    st->lights.type.push_back(ygl::gl_light_type::point);
                }
            }
        }
    } else {
        for (auto ist : scn->meshes) {
            for (auto shp : ist->shapes) {
                if (!shp->mat) continue;
                auto mat = shp->mat;
                if (mat->emission == ygl::zero3f) continue;
                for (auto p : shp->points) {
                    if (st->lights.pos.size() >= 16) break;
                    st->lights.pos.push_back(shp->pos[p]);
                    st->lights.ke.push_back(mat->emission);
                    st->lights.type.push_back(ygl::gl_light_type::point);
                }
            }
        }
    }
}

// Init shading
inline void update_shade_state(
    const ygl::gltf_scene_group* sc, shade_state* st) {
    if (!is_program_valid(st->prog)) st->prog = ygl::make_stdsurface_program();
    st->txt[nullptr] = {};
    for (auto txt : sc->textures) {
        if (st->txt.find(txt) != st->txt.end()) continue;
        if (!txt->hdr.empty()) {
            st->txt[txt] = ygl::make_texture(txt->hdr, true, true, true);
        } else if (!txt->ldr.empty()) {
            st->txt[txt] = ygl::make_texture(txt->ldr, true, true, true);
        } else
            assert(false);
    }
    for (auto mesh : sc->meshes) {
        for (auto shape : mesh->shapes) {
            if (st->vbo.find(shape) != st->vbo.end()) continue;
            st->vbo[shape] = shape_vbo();
            if (!shape->pos.empty())
                st->vbo[shape].pos = ygl::make_vertex_buffer(shape->pos);
            if (!shape->norm.empty())
                st->vbo[shape].norm = ygl::make_vertex_buffer(shape->norm);
            if (!shape->texcoord.empty())
                st->vbo[shape].texcoord =
                    ygl::make_vertex_buffer(shape->texcoord);
            if (!shape->texcoord1.empty())
                st->vbo[shape].texcoord1 =
                    ygl::make_vertex_buffer(shape->texcoord1);
            if (!shape->color.empty())
                st->vbo[shape].color = ygl::make_vertex_buffer(shape->color);
            if (!shape->tangsp.empty())
                st->vbo[shape].tangsp = ygl::make_vertex_buffer(shape->tangsp);
            if (!shape->skin_weights.empty())
                st->vbo[shape].skin_weights =
                    ygl::make_vertex_buffer(shape->skin_weights);
            if (!shape->skin_joints.empty())
                st->vbo[shape].skin_joints =
                    ygl::make_vertex_buffer(shape->skin_joints);
            if (!shape->points.empty())
                st->vbo[shape].points = ygl::make_element_buffer(shape->points);
            if (!shape->lines.empty())
                st->vbo[shape].lines = ygl::make_element_buffer(shape->lines);
            if (!shape->triangles.empty())
                st->vbo[shape].triangles =
                    ygl::make_element_buffer(shape->triangles);
        }
    }
}

// Draw a mesh
inline void shade_mesh(const ygl::gltf_mesh* msh, const ygl::gltf_skin* sk,
    const std::vector<float>& morph_weights, shade_state* st,
    const ygl::mat4f& xform, bool highlighted, bool edges, bool wireframe,
    bool cutout) {
    static auto default_material = ygl::gltf_material();
    default_material.metallic_roughness =
        new ygl::gltf_material_metallic_roughness();
    default_material.metallic_roughness->base = {0.2f, 0.2f, 0.2f};

    auto txt_info =
        [st](ygl::gltf_texture* gtxt,
            const ygl::gltf_texture_info* ginfo) -> ygl::gl_texture_info {
        auto info = ygl::gl_texture_info();
        if (!gtxt) return info;
        info.txt = st->txt.at(gtxt);
        if (!ginfo) return info;
        info.wrap_s = (ygl::gl_texture_wrap)ginfo->wrap_s;
        info.wrap_t = (ygl::gl_texture_wrap)ginfo->wrap_t;
        info.filter_min = (ygl::gl_texture_filter)ginfo->filter_min;
        info.filter_mag = (ygl::gl_texture_filter)ginfo->filter_mag;
        return info;
    };

    for (auto shp : msh->shapes) {
        begin_stdsurface_shape(st->prog, xform);

        auto etype = ygl::gl_elem_type::triangle;
        if (!shp->lines.empty()) etype = ygl::gl_elem_type::line;
        if (!shp->points.empty()) etype = ygl::gl_elem_type::point;

        set_stdsurface_highlight(
            st->prog, (highlighted) ? ygl::vec4f{1, 1, 0, 1} : ygl::zero4f);

        auto mat = (shp->mat) ? shp->mat : &default_material;
        float op = 1;
        if (mat->specular_glossiness) {
            auto sg = mat->specular_glossiness;
            op = sg->opacity;
            set_stdsurface_material(st->prog,
                ygl::material_type::specular_glossiness, etype, mat->emission,
                sg->diffuse, sg->specular, sg->glossiness, sg->opacity,
                txt_info(mat->emission_txt, mat->emission_txt_info),
                txt_info(sg->diffuse_txt, sg->diffuse_txt_info),
                txt_info(sg->specular_txt, sg->specular_txt_info), {},
                txt_info(mat->normal_txt, mat->normal_txt_info),
                txt_info(mat->occlusion_txt, mat->occlusion_txt_info), false,
                mat->double_sided, cutout);
        } else if (mat->metallic_roughness) {
            auto mr = mat->metallic_roughness;
            op = mr->opacity;
            set_stdsurface_material(st->prog,
                ygl::material_type::metallic_roughness, etype, mat->emission,
                mr->base, {mr->metallic, mr->metallic, mr->metallic},
                mr->roughness, mr->opacity,
                txt_info(mat->emission_txt, mat->emission_txt_info),
                txt_info(mr->base_txt, mr->base_txt_info),
                txt_info(mr->metallic_txt, mr->metallic_txt_info), {},
                txt_info(mat->normal_txt, mat->normal_txt_info),
                txt_info(mat->occlusion_txt, mat->occlusion_txt_info), false,
                mat->double_sided, cutout);
        } else {
            set_stdsurface_material(st->prog,
                ygl::material_type::specular_roughness, etype, mat->emission,
                ygl::zero3f, ygl::zero3f, 0.5f, 1,
                txt_info(mat->emission_txt, mat->emission_txt_info), {}, {}, {},
                {}, {}, false, mat->double_sided, cutout);
        }

        auto vbo = st->vbo.at(shp);
        set_stdsurface_vert(
            st->prog, vbo.pos, vbo.norm, vbo.texcoord, vbo.color, vbo.tangsp);
        if (sk) {
            auto skin_xforms = ygl::get_skin_transforms(sk, xform);
#if 0
            st->prog.set_vert_gltf_skinning(st->prog,
                                                    vert_vbo.skin_weights, vert_vbo.skin_joints, skin_xforms.size(),
                                                    skin_xforms.data());
#else
            std::vector<ygl::vec3f> skinned_pos, skinned_norm;
            // ygl::compute_skinning(shp->pos, shp->norm, shp->skin_weights,
            //     shp->skin_joints, skin_xforms, skinned_pos, skinned_norm);
            ygl::compute_matrix_skinning(shp->pos, shp->norm, shp->skin_weights,
                shp->skin_joints, skin_xforms, skinned_pos, skinned_norm);
            update_vertex_buffer(vbo.pos, skinned_pos);
            update_vertex_buffer(vbo.norm, skinned_norm);
#endif
        } else if (!morph_weights.empty() && !shp->morph_targets.empty()) {
            std::vector<ygl::vec3f> morph_pos, morph_norm;
            std::vector<ygl::vec4f> morph_tang;
            ygl::compute_morphing_deformation(
                shp, morph_weights, morph_pos, morph_norm, morph_tang);
            update_vertex_buffer(vbo.pos, morph_pos);
            update_vertex_buffer(vbo.norm, morph_norm);
        } else {
            set_stdsurface_vert_skinning_off(st->prog);
        }

        ygl::gl_enable_culling(!mat->double_sided);

        draw_elems(vbo.points);
        draw_elems(vbo.lines);
        draw_elems(vbo.triangles);

        ygl::gl_enable_culling(false);

        if (edges && !wireframe) {
            assert(ygl::gl_check_error());
            set_stdsurface_material(st->prog,
                ygl::material_type::specular_roughness, etype, ygl::zero3f,
                ygl::zero3f, ygl::zero3f, 0.5f, 1, {}, {}, {}, {}, {}, {}, true,
                mat->double_sided, cutout);

            assert(ygl::gl_check_error());
            ygl::gl_line_width(2);
            draw_elems(vbo.points);
            draw_elems(vbo.lines);
            draw_elems(vbo.triangles);
            ygl::gl_line_width(1);
            assert(ygl::gl_check_error());
        }

        end_stdsurface_shape(st->prog);
    }
}

// Display a scene
inline void shade_scene(const ygl::gltf_scene_group* scns, shade_state* st,
    const ygl::camera* ycam, const ygl::gltf_node* gcam, void* selection,
    const ygl::vec4f& background, float exposure, float gamma, bool filmic,
    bool wireframe, bool edges, bool cutout, bool camera_lights,
    const ygl::vec3f& amb) {
    // update state
    update_shade_state(scns, st);

    // begin frame
    ygl::gl_enable_depth_test(true);
    ygl::gl_enable_culling(false);

    ygl::gl_enable_wireframe(wireframe);

    ygl::mat4f camera_xform, camera_view, camera_proj;
    if (gcam) {
        camera_xform = ygl::mat4f(gcam->xform());
        camera_view = ygl::inverse(ygl::mat4f(gcam->xform()));
        if (gcam->cam->ortho) {
            auto near = (gcam->cam->near) ? gcam->cam->near : 0.001f;
            auto far = (gcam->cam->far) ? gcam->cam->far : 10000;
            camera_proj = ygl::ortho2d_mat(gcam->cam->yfov * gcam->cam->aspect,
                gcam->cam->yfov, near, far);
        } else {
            auto near = (gcam->cam->near) ? gcam->cam->near : 0.001f;
            if (gcam->cam->far) {
                camera_proj = ygl::perspective_mat(
                    gcam->cam->yfov, gcam->cam->aspect, near, gcam->cam->far);
            } else {
                camera_proj = ygl::perspective_mat(
                    gcam->cam->yfov, gcam->cam->aspect, 0.01f);
            }
        }
    } else {
        camera_xform = ygl::frame_to_mat(ycam->frame);
        camera_view = ygl::frame_to_mat(ygl::inverse(ycam->frame));
        auto near = (ycam->near) ? ycam->near : 0.001f;
        if (ycam->far) {
            camera_proj =
                ygl::perspective_mat(ycam->yfov, ycam->aspect, near, ycam->far);
        } else {
            camera_proj = ygl::perspective_mat(ycam->yfov, ycam->aspect, near);
        }
    }

    begin_stdsurface_frame(st->prog, camera_lights, exposure, gamma, filmic,
        camera_xform, camera_view, camera_proj);

    if (!camera_lights) {
        update_shade_lights(st, scns);
        set_stdsurface_lights(st->prog, amb, st->lights);
    }

    auto instances = ygl::get_mesh_nodes(scns->default_scene);

    if (!instances.empty()) {
        for (auto ist : instances) {
            shade_mesh(ist->msh, ist->skn, ist->morph_weights, st, ist->xform(),
                ist == selection || ist->msh == selection, edges, wireframe,
                cutout);
        }
    } else {
        for (auto msh : scns->meshes) {
            shade_mesh(msh, nullptr, {}, st, ygl::identity_mat4f,
                msh == selection, edges, wireframe, cutout);
        }
    }

    end_stdsurface_frame(st->prog);

    ygl::gl_enable_wireframe(false);
}

void draw_scene(ygl::gl_window* win) {
    auto scn = (app_state*)get_user_pointer(win);
    auto window_size = get_window_size(win);
    auto aspect = (float)window_size.x / (float)window_size.y;
    scn->view_cam->aspect = aspect;
    if (scn->gcam) scn->gcam->cam->aspect = aspect;
    shade_scene(scn->gscn, scn->shstate, scn->view_cam, scn->gcam,
        scn->selection, scn->background, scn->exposure, scn->gamma, scn->filmic,
        scn->wireframe, scn->edges, scn->alpha_cutout, scn->camera_lights,
        scn->amb);
}

void draw_tree_widgets(ygl::gl_window* win, const std::string& lbl,
    ygl::gltf_camera* cam, void*& selection) {
    draw_tree_widget_leaf(win, lbl + cam->name, selection, cam);
}

void draw_tree_widgets(ygl::gl_window* win, const std::string& lbl,
    ygl::gltf_texture* txt, void*& selection) {
    draw_tree_widget_leaf(win, lbl + txt->path, selection, txt);
}

void draw_tree_widgets(ygl::gl_window* win, const std::string& lbl,
    ygl::gltf_material_metallic_roughness* mat, void*& selection) {
    if (mat->base_txt) draw_tree_widgets(win, "kb: ", mat->base_txt, selection);
    if (mat->metallic_txt)
        draw_tree_widgets(win, "km: ", mat->metallic_txt, selection);
}

void draw_tree_widgets(ygl::gl_window* win, const std::string& lbl,
    ygl::gltf_material_specular_glossiness* mat, void*& selection) {
    if (mat->diffuse_txt)
        draw_tree_widgets(win, "kd: ", mat->diffuse_txt, selection);
    if (mat->specular_txt)
        draw_tree_widgets(win, "ks: ", mat->specular_txt, selection);
}

void draw_tree_widgets(ygl::gl_window* win, const std::string& lbl,
    ygl::gltf_material* mat, void*& selection) {
    if (draw_tree_widget_begin(win, lbl + mat->name, selection, mat)) {
        if (mat->emission_txt)
            draw_tree_widgets(win, "ke: ", mat->emission_txt, selection);
        if (mat->metallic_roughness)
            draw_tree_widgets(win, "mr: ", mat->metallic_roughness, selection);
        if (mat->specular_glossiness)
            draw_tree_widgets(win, "sg: ", mat->specular_glossiness, selection);
        if (mat->normal_txt)
            draw_tree_widgets(win, "norm: ", mat->normal_txt, selection);
        if (mat->occlusion_txt)
            draw_tree_widgets(win, "occ: ", mat->occlusion_txt, selection);
        draw_tree_widget_end(win);
    }
}

void draw_tree_widgets(ygl::gl_window* win, const std::string& lbl,
    ygl::gltf_shape* shp, void*& selection) {
    if (draw_tree_widget_begin(win, lbl + shp->name, selection, shp)) {
        if (shp->mat) draw_tree_widgets(win, "mat: ", shp->mat, selection);
        draw_tree_widget_end(win);
    }
}

void draw_tree_widgets(ygl::gl_window* win, const std::string& lbl,
    ygl::gltf_mesh* msh, void*& selection) {
    if (draw_tree_widget_begin(win, lbl + msh->name, selection, msh)) {
        auto sid = 0;
        for (auto shp : msh->shapes) {
            draw_tree_widgets(
                win, "sh" + std::to_string(sid++) + ": ", shp, selection);
        }
        draw_tree_widget_end(win);
    }
}

void draw_tree_widgets(ygl::gl_window* win, const std::string& lbl,
    ygl::gltf_node* node, void*& selection) {
    if (draw_tree_widget_begin(win, lbl + node->name, selection, node)) {
        if (node->msh) draw_tree_widgets(win, "mesh: ", node->msh, selection);
        if (node->cam) draw_tree_widgets(win, "cam: ", node->cam, selection);
        for (auto child : node->children)
            draw_tree_widgets(win, "", child, selection);
        draw_tree_widget_end(win);
    }
}

void draw_tree_widgets(ygl::gl_window* win, const std::string& lbl,
    ygl::gltf_scene* scn, void*& selection) {
    if (draw_tree_widget_begin(win, lbl + scn->name, selection, scn)) {
        for (auto node : scn->nodes)
            draw_tree_widgets(win, "", node, selection);
        draw_tree_widget_end(win);
    }
}

void draw_tree_widgets(ygl::gl_window* win, const std::string& lbl,
    ygl::gltf_animation* anim, void*& selection) {
    if (draw_tree_widget_begin(win, lbl, selection, anim)) {
        auto sid = 0;
        for (auto node : anim->nodes) {
            draw_tree_widgets(
                win, "node " + std::to_string(sid++) + ": ", node, selection);
        }
        draw_tree_widget_end(win);
    }
}

void draw_tree_widgets(ygl::gl_window* win, const std::string& lbl,
    ygl::gltf_animation_group* anims, void*& selection) {
    if (draw_tree_widget_begin(win, lbl + anims->name, selection, anims)) {
        auto sid = 0;
        for (auto anim : anims->animations) {
            draw_tree_widgets(
                win, "anim " + std::to_string(sid++), anim, selection);
        }
        draw_tree_widget_end(win);
    }
}

void draw_tree_widgets(ygl::gl_window* win, const std::string& lbl,
    ygl::gltf_skin* skin, void*& selection) {
    if (draw_tree_widget_begin(win, lbl + skin->name, selection, skin)) {
        draw_tree_widgets(win, "skeleton", skin->root, selection);
        draw_tree_widget_end(win);
    }
}

void draw_tree_widgets(ygl::gl_window* win, const std::string& lbl,
    ygl::gltf_scene_group* oscn, void*& selection) {
    if (draw_tree_widget_begin(win, lbl + "scenes")) {
        for (auto scn : oscn->scenes)
            draw_tree_widgets(win, "", scn, selection);
        draw_tree_widget_end(win);
    }

    if (draw_tree_widget_begin(win, lbl + "cameras")) {
        for (auto cam : oscn->cameras)
            draw_tree_widgets(win, "", cam, selection);
        draw_tree_widget_end(win);
    }

    if (draw_tree_widget_begin(win, lbl + "meshes")) {
        for (auto msh : oscn->meshes)
            draw_tree_widgets(win, "", msh, selection);
        draw_tree_widget_end(win);
    }

    if (draw_tree_widget_begin(win, lbl + "nodes")) {
        for (auto node : oscn->nodes) {
            if (!node->parent) draw_tree_widgets(win, "", node, selection);
        }
        draw_tree_widget_end(win);
    }

    if (draw_tree_widget_begin(win, lbl + "materials")) {
        for (auto mat : oscn->materials)
            draw_tree_widgets(win, "", mat, selection);
        draw_tree_widget_end(win);
    }

    if (draw_tree_widget_begin(win, lbl + "textures")) {
        for (auto txt : oscn->textures)
            draw_tree_widgets(win, "", txt, selection);
        draw_tree_widget_end(win);
    }

    if (draw_tree_widget_begin(win, lbl + "animations")) {
        for (auto anim : oscn->animations)
            draw_tree_widgets(win, "", anim, selection);
        draw_tree_widget_end(win);
    }

    if (draw_tree_widget_begin(win, lbl + "skins")) {
        for (auto skin : oscn->skins)
            draw_tree_widgets(win, "", skin, selection);
        draw_tree_widget_end(win);
    }
}

void draw_elem_widgets(ygl::gl_window* win, ygl::gltf_scene_group* gscn,
    ygl::gltf_texture* txt, void*& selection, const shade_state* state) {
    if (selection != txt) return;

    draw_separator_widget(win);
    draw_label_widget(win, "path", txt->path);
    auto str = ygl::format("{} x {} @ 4 {}", txt->width(), txt->height(),
        (!txt->ldr.empty()) ? "byte" : "float");
    draw_label_widget(win, "size", str);
    if (state->txt.find(txt) != state->txt.end()) {
        draw_image_widget(win, get_texture_id(state->txt.at(txt)), {128, 128},
            {txt->width(), txt->height()});
    }
}

void draw_elem_widgets(ygl::gl_window* win, ygl::gltf_scene_group* gscn,
    ygl::gltf_texture_info* info, void*& selection, const shade_state* state) {
    static const auto wrap_names =
        std::vector<std::pair<std::string, ygl::gltf_texture_wrap>>{
            {"repeat", ygl::gltf_texture_wrap::repeat},
            {"clamp", ygl::gltf_texture_wrap::clamp},
            {"mirror", ygl::gltf_texture_wrap::mirror}};
    static const auto filter_mag_names =
        std::vector<std::pair<std::string, ygl::gltf_texture_filter>>{
            {"linear", ygl::gltf_texture_filter::linear},
            {"nearest", ygl::gltf_texture_filter::nearest},
        };
    static const auto filter_min_names =
        std::vector<std::pair<std::string, ygl::gltf_texture_filter>>{
            {"linear", ygl::gltf_texture_filter::linear},
            {"nearest", ygl::gltf_texture_filter::nearest},
            {"linear_mipmap_linear",
                ygl::gltf_texture_filter::linear_mipmap_linear},
            {"linear_mipmap_nearest",
                ygl::gltf_texture_filter::linear_mipmap_nearest},
            {"nearest_mipmap_linear",
                ygl::gltf_texture_filter::nearest_mipmap_linear},
            {"nearest_mipmap_nearest",
                ygl::gltf_texture_filter::nearest_mipmap_nearest},
        };

    if (selection != info) return;

    draw_indent_widget_begin(win);

    draw_combo_widget(win, "wrap s", info->wrap_s, wrap_names);
    draw_combo_widget(win, "wrap t", info->wrap_t, wrap_names);
    draw_combo_widget(win, "filter mag", info->filter_mag, filter_mag_names);
    draw_combo_widget(win, "filter min", info->filter_min, filter_min_names);

    draw_indent_widget_end(win);
}

void draw_elem_widgets(ygl::gl_window* win, ygl::gltf_scene_group* gscn,
    ygl::gltf_material_metallic_roughness* mat, void*& selection,
    const shade_state* state) {
    auto txt_names = std::vector<std::pair<std::string, ygl::gltf_texture*>>{
        {"<none>", nullptr}};
    for (auto txt : gscn->textures) txt_names.push_back({txt->path, txt});

    draw_value_widget(win, "mr base", mat->base, 0, 1);
    draw_value_widget(win, "mr opacity", mat->opacity, 0, 1);
    draw_value_widget(win, "mr metallic", mat->metallic, 0, 1);
    draw_value_widget(win, "mr roughness", mat->roughness, 0, 1);
    draw_combo_widget(win, "mr base txt", mat->base_txt, txt_names);
    if (mat->base_txt && mat->base_txt_info)
        draw_elem_widgets(win, gscn, mat->base_txt_info, selection, state);
    draw_combo_widget(win, "mr metallic txt", mat->metallic_txt, txt_names);
    if (mat->metallic_txt && mat->metallic_txt_info)
        draw_elem_widgets(win, gscn, mat->metallic_txt_info, selection, state);
}

void draw_elem_widgets(ygl::gl_window* win, ygl::gltf_scene_group* gscn,
    ygl::gltf_material_specular_glossiness* mat, void*& selection,
    const shade_state* state) {
    auto txt_names = std::vector<std::pair<std::string, ygl::gltf_texture*>>{
        {"<none>", nullptr}};
    for (auto txt : gscn->textures) txt_names.push_back({txt->path, txt});

    draw_value_widget(win, "sg diffuse", mat->diffuse, 0, 1);
    draw_value_widget(win, "sg opacity", mat->opacity, 0, 1);
    draw_value_widget(win, "sg specular", mat->specular, 0, 1);
    draw_value_widget(win, "sg glossiness", mat->glossiness, 0, 1);
    draw_combo_widget(win, "sg diffuse txt", mat->diffuse_txt, txt_names);
    if (mat->diffuse_txt && mat->diffuse_txt_info)
        draw_elem_widgets(win, gscn, mat->diffuse_txt_info, selection, state);
    draw_combo_widget(win, "sg specular txt", mat->specular_txt, txt_names);
    if (mat->specular_txt && mat->specular_txt_info)
        draw_elem_widgets(win, gscn, mat->specular_txt_info, selection, state);
}

void draw_elem_widgets(ygl::gl_window* win, ygl::gltf_scene_group* gscn,
    ygl::gltf_material* mat, void*& selection, const shade_state* state) {
    if (selection != mat) return;

    auto txt_names = std::vector<std::pair<std::string, ygl::gltf_texture*>>{
        {"<none>", nullptr}};
    for (auto txt : gscn->textures) txt_names.push_back({txt->path, txt});

    draw_separator_widget(win);
    draw_label_widget(win, "name", mat->name);
    draw_value_widget(win, "emission", mat->emission, 0, 1000);
    draw_combo_widget(win, "emission txt", mat->emission_txt, txt_names);
    if (mat->emission_txt && mat->emission_txt_info)
        draw_elem_widgets(win, gscn, mat->emission_txt_info, selection, state);
    draw_combo_widget(win, "normal txt", mat->normal_txt, txt_names);
    if (mat->normal_txt && mat->normal_txt_info)
        draw_elem_widgets(win, gscn, mat->normal_txt_info, selection, state);
    draw_combo_widget(win, "occlusion txt", mat->occlusion_txt, txt_names);
    if (mat->occlusion_txt && mat->occlusion_txt_info)
        draw_elem_widgets(win, gscn, mat->occlusion_txt_info, selection, state);
    if (mat->metallic_roughness) {
        draw_separator_widget(win);
        draw_elem_widgets(win, gscn, mat->metallic_roughness, selection, state);
    } else {
        if (draw_button_widget(win, "add metallic roughness")) {
            mat->metallic_roughness =
                new ygl::gltf_material_metallic_roughness();
        }
    }
    if (mat->specular_glossiness) {
        draw_separator_widget(win);
        draw_elem_widgets(
            win, gscn, mat->specular_glossiness, selection, state);
    } else {
        if (draw_button_widget(win, "add specular glossiness")) {
            mat->specular_glossiness =
                new ygl::gltf_material_specular_glossiness();
        }
    }
}

void draw_elem_widgets(ygl::gl_window* win, ygl::gltf_scene_group* gscn,
    ygl::gltf_shape* shp, void*& selection, const shade_state* state) {
    if (selection != shp) return;

    auto mat_names = std::vector<std::pair<std::string, ygl::gltf_material*>>{
        {"<none>", nullptr}};
    for (auto mat : gscn->materials) mat_names.push_back({mat->name, mat});

    draw_separator_widget(win);
    draw_label_widget(win, "name", shp->name);
    draw_combo_widget(win, "material", shp->mat, mat_names);
    draw_label_widget(win, "verts", (int)shp->pos.size());
    if (!shp->triangles.empty())
        draw_label_widget(win, "triangles", (int)shp->triangles.size());
    if (!shp->lines.empty())
        draw_label_widget(win, "lines", (int)shp->lines.size());
    if (!shp->points.empty())
        draw_label_widget(win, "points", (int)shp->points.size());
}

void draw_elem_widgets(ygl::gl_window* win, ygl::gltf_scene_group* gscn,
    ygl::gltf_mesh* msh, void*& selection, const shade_state* state) {
    if (selection != msh) return;
    draw_separator_widget(win);
    draw_label_widget(win, "name", msh->name);
    for (auto shp : msh->shapes) {
        draw_elem_widgets(win, gscn, shp, selection, state);
    }
}

void draw_elem_widgets(ygl::gl_window* win, ygl::gltf_scene_group* gscn,
    ygl::gltf_camera* cam, void*& selection, const shade_state* state) {
    if (selection != cam) return;
    draw_separator_widget(win);
    draw_label_widget(win, "name", cam->name);
    draw_value_widget(win, "yfov", cam->yfov, 0.1, 4);
    draw_value_widget(win, "aspect", cam->aspect, 0.1, 4);
    draw_value_widget(win, "focus", cam->focus, 0.01, 10);
    draw_value_widget(win, "aperture", cam->aperture, 0, 1);
}

void draw_elem_widgets(ygl::gl_window* win, ygl::gltf_scene_group* gscn,
    ygl::gltf_animation* anim, void*& selection, const shade_state* state) {
    if (selection != anim) return;
    auto interp_name = std::map<ygl::gltf_animation_interpolation, std::string>{
        {ygl::gltf_animation_interpolation::step, "step"},
        {ygl::gltf_animation_interpolation::linear, "linear"},
        {ygl::gltf_animation_interpolation::catmull_rom, "catmull-rom spline"},
        {ygl::gltf_animation_interpolation::cubic, "cubic spline"},
    };
    auto node_names = std::vector<std::pair<std::string, ygl::gltf_node*>>{
        {"<none>", nullptr}};
    for (auto node : gscn->nodes) node_names.push_back({node->name, node});
    draw_separator_widget(win);
    draw_label_widget(win, "interpolation", interp_name.at(anim->interp));
    draw_label_widget(win, "num times", (int)anim->time.size());
    draw_label_widget(win, "num translation", (int)anim->translation.size());
    draw_label_widget(win, "num rotation", (int)anim->rotation.size());
    draw_label_widget(win, "num scale", (int)anim->scale.size());
    auto selected_node = (ygl::gltf_node*)nullptr;
    draw_combo_widget(win, "", selected_node, node_names);
}

void draw_elem_widgets(ygl::gl_window* win, ygl::gltf_scene_group* gscn,
    ygl::gltf_animation_group* anims, void*& selection,
    const shade_state* state) {
    if (selection != anims) return;
    draw_separator_widget(win);
    draw_label_widget(win, "name", anims->name);
    for (auto anim : anims->animations) {
        draw_elem_widgets(win, gscn, anim, selection, state);
    }
}

void draw_elem_widgets(ygl::gl_window* win, ygl::gltf_scene_group* gscn,
    ygl::gltf_skin* skin, void*& selection, const shade_state* state) {
    if (selection != skin) return;
    draw_separator_widget(win);
    draw_label_widget(win, "name", skin->name);
}

void draw_elem_widgets(ygl::gl_window* win, ygl::gltf_scene_group* gscn,
    ygl::gltf_node* node, void*& selection, const shade_state* state) {
    if (selection != node) return;

    auto msh_names = std::vector<std::pair<std::string, ygl::gltf_mesh*>>{
        {"<none>", nullptr}};
    for (auto msh : gscn->meshes) msh_names.push_back({msh->name, msh});
    auto cam_names = std::vector<std::pair<std::string, ygl::gltf_camera*>>{
        {"<none>", nullptr}};
    for (auto cam : gscn->cameras) cam_names.push_back({cam->name, cam});

    draw_separator_widget(win);
    draw_label_widget(win, "name", node->name);
    draw_combo_widget(win, "mesh", node->msh, msh_names);
    draw_combo_widget(win, "cam", node->cam, cam_names);
    draw_value_widget(win, "translation", node->translation, -10, 10);
    draw_value_widget(win, "rotation", node->rotation);
    draw_value_widget(win, "scale", node->scale, 0.01f, 10);
    draw_value_widget(win, "matrix", node->matrix, -10, 10);
}

void draw_elem_widgets(ygl::gl_window* win, ygl::gltf_scene_group* gscn,
    void*& selection, const shade_state* state) {
    for (auto cam : gscn->cameras) {
        draw_elem_widgets(win, gscn, cam, selection, state);
    }

    for (auto msh : gscn->meshes) {
        draw_elem_widgets(win, gscn, msh, selection, state);
        for (auto shp : msh->shapes) {
            draw_elem_widgets(win, gscn, shp, selection, state);
        }
    }

    for (auto node : gscn->nodes) {
        draw_elem_widgets(win, gscn, node, selection, state);
    }

    for (auto mat : gscn->materials) {
        draw_elem_widgets(win, gscn, mat, selection, state);
    }

    for (auto txt : gscn->textures) {
        draw_elem_widgets(win, gscn, txt, selection, state);
    }

    for (auto anims : gscn->animations) {
        draw_elem_widgets(win, gscn, anims, selection, state);
        for (auto anim : anims->animations) {
            draw_elem_widgets(win, gscn, anim, selection, state);
        }
    }

    for (auto skin : gscn->skins) {
        draw_elem_widgets(win, gscn, skin, selection, state);
    }
}

void draw_scene_widgets(ygl::gl_window* win, ygl::gltf_scene_group* gscn,
    void*& selection, const shade_state* state) {
    draw_scroll_widget_begin(win, "model", 240, false);
    draw_tree_widgets(win, "", gscn, selection);
    draw_scroll_widget_end(win);
    draw_elem_widgets(win, gscn, selection, state);
}

#ifndef YSCENE_SIMULATE_HACK
void simulate_step(app_state* scene) {}
#else
void simulate_step(app_state* scene);
#endif

void draw_widgets(ygl::gl_window* win) {
    auto scn = (app_state*)get_user_pointer(win);
    if (begin_widgets(win, "yview")) {
        draw_label_widget(win, "scene", scn->filename);
        if (scn->time_range != ygl::zero2f) {
            draw_value_widget(
                win, "time", scn->time, scn->time_range.x, scn->time_range.y);
            draw_value_widget(win, "play animation", scn->animate);
        }
        if (scn->gscn->default_scene) {
            auto camera_names =
                std::vector<std::pair<std::string, ygl::gltf_node*>>{
                    {"<view>", nullptr}};
            for (auto cam : ygl::get_camera_nodes(scn->gscn->default_scene)) {
                camera_names.push_back({cam->name, cam});
            }
            if (draw_combo_widget(win, "camera", scn->gcam, camera_names))
                scn->scene_updated = true;
        }
        draw_value_widget(win, "wire", scn->wireframe);
        draw_value_widget(win, "edges", scn->edges);
        draw_value_widget(win, "cutout", scn->alpha_cutout);
        draw_value_widget(win, "fps", scn->navigation_fps);
        draw_value_widget(win, "hdr exposure", scn->exposure, -20, 20, 1);
        draw_value_widget(win, "hdr gamma", scn->gamma, 0.1, 5);
        draw_value_widget(win, "hdr filmic", scn->filmic);
        if (draw_header_widget(win, "view cam")) {
            auto cam = scn->view_cam;
            draw_value_widget(win, "yfov", cam->yfov, 0.1, 5);
            draw_value_widget(win, "aspect", cam->aspect, 0.1, 3);
            draw_value_widget(win, "focus", cam->focus, 0.1, 1000);
            draw_value_widget(win, "near", cam->near, 0.01, 1);
            draw_value_widget(win, "far", cam->far, 1, 10000);
        }
        if (draw_header_widget(win, "inspect")) {
            draw_scene_widgets(win, scn->gscn, scn->selection, scn->shstate);
        }
    }
    end_widgets(win);
}

void parse_cmdline(app_state* scene, int argc, char** argv, const char* name,
    const char* help, bool simulation, bool trace) {
    using namespace std::string_literals;

    // parse command line
    auto parser = ygl::make_parser(argc, argv, name, help);

    // render
    scene->exposure =
        parse_opt(parser, "--exposure", "-e", "hdr image exposure", 0.0f);
    scene->gamma = parse_opt(parser, "--gamma", "-g", "hdr image gamma", 2.2f);
    scene->filmic = parse_flag(parser, "--filmic", "-F", "hdr filmic output");
    scene->resolution =
        parse_opt(parser, "--resolution", "-r", "image resolution", 540);

    auto amb = parse_opt(parser, "--ambient", "", "ambient factor", 0.0f);

    scene->camera_lights = parse_flag(
        parser, "--camera-lights", "-c", "enable camera lights", false);

    scene->amb = {amb, amb, amb};

    // load options
    scene->split_shapes = parse_flag(
        parser, "--split-shapes", "", "split meshes into single shapes", false);
    scene->add_spec_gloss = parse_flag(
        parser, "--add-spec-gloss", "", "add spec-gloss to gltf", false);

    // params
    scene->imfilename =
        parse_opt(parser, "--output-image", "-o", "image filename", "out.hdr"s);
    scene->filename = parse_arg(parser, "scene", "scene filename", ""s);

    // check parsing
    if (should_exit(parser)) {
        printf("%s\n", get_usage(parser).c_str());
        exit(1);
    }
}

// init draw with shading
void shade_init(ygl::gl_window* win) {
    auto scn = (app_state*)get_user_pointer(win);
    scn->shstate = new shade_state();
    update_shade_state(scn->gscn, scn->shstate);
}

// draw with shading
void shade_draw(ygl::gl_window* win) {
    auto scn = (app_state*)get_user_pointer(win);
    if (scn->gscn) {
        ygl::update_animated_transforms(scn->gscn, scn->time);
        ygl::update_transforms(scn->gscn);
    }

    auto fwh = get_framebuffer_size(win);
    ygl::gl_set_viewport(fwh);

    ygl::gl_clear_buffers();
    draw_scene(win);
    draw_widgets(win);
    swap_buffers(win);
    if (scn->shstate->lights.pos.empty()) scn->camera_lights = true;
}

bool update(app_state* scn) {
    // advance time
    if (scn->animate && scn->time_range != ygl::zero2f) {
        if (scn->time >= scn->time_range.y) {
            scn->time = scn->time_range.x;
        } else {
            scn->time += 1 / 60.0f;
        }
        return true;
    }
    return false;
}

// run ui loop
void run_ui(app_state* app, int w, int h) {
    // window
    auto win = ygl::make_window(w, h, "ygltfview | " + app->filename, app);
    set_window_callbacks(win, nullptr, nullptr, shade_draw);

    // load textures
    shade_init(win);

    // init widget
    init_widgets(win);

    // loop
    while (!should_close(win)) {
        // handle mouse and keyboard for navigation
        if (!app->gcam) {
            handle_camera_navigation(win, app->view_cam, app->navigation_fps);
        }

        // draw
        shade_draw(win);

        // update
        update(app);

        // event hadling
        poll_events(win);
    }

    clear_window(win);
}

void draw_custom_widgets(ygl::gl_window* win) {
    auto scn = (app_state*)get_user_pointer(win);
    if (scn->time_range != ygl::zero2f) {
        draw_value_widget(
            win, "time", scn->time, scn->time_range.x, scn->time_range.y);
        draw_value_widget(win, "play animation", scn->animate);
    }
}

int main(int argc, char* argv[]) {
    // create empty scene
    auto scn = new app_state();

    // params
    parse_cmdline(
        scn, argc, argv, "yview", "interactively view scenes", false, false);

    // scene
    log_info("loading scene {}", scn->filename);
    if (!load_scene(scn, scn->filename, false, false)) return 1;

    // run ui
    auto width = (int)std::round(scn->view_cam->aspect * scn->resolution);
    auto height = scn->resolution;
    run_ui(scn, width, height);

    // clear
    delete scn;

    // done
    return 0;
}
