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

#ifndef _YSCENE_H_

#include "../yocto/yocto_gltf.h"
#include "../yocto/yocto_img.h"
#include "../yocto/yocto_math.h"
#include "../yocto/yocto_obj.h"
#include "../yocto/yocto_sym.h"
#include "../yocto/yocto_trace.h"
#include "../yocto/yocto_utils.h"

#ifndef YOCTO_NO_OPENGL
#include "../yocto/yocto_glu.h"
#include "../yocto/yocto_gui.h"
#endif

// ---------------------------------------------------------------------------
// SCENE (OBJ or GLTF) AND APPLICATION PARAMETERS
// ---------------------------------------------------------------------------

#ifndef YOCTO_NO_OPENGL

// convenience typedef
using uint = unsigned int;

//
// OpenGL state
//
struct yshade_state {
    struct vert_vbo {
        uint pos, norm, texcoord, texcoord1, color, tangsp, skin_joints,
            skin_weights;
    };
    struct elem_vbo {
        uint points, lines, triangles;
    };

    // shade state
    uint prog = 0;
    uint vao = 0;
    std::map<void*, uint> txt;
    std::map<void*, vert_vbo> vert;
    std::map<void*, elem_vbo> elem;

    // lights
    std::vector<ym::vec3f> lights_pos;
    std::vector<ym::vec3f> lights_ke;
    std::vector<yglu::ltype> lights_ltype;
};
#else
struct yshade_state {};
#endif

//
// Camera
//
struct ycamera {
    ym::frame3f frame = ym::identity_frame3f;
    float yfov = 1;
    float aspect = 1;
    float focus = 1;
    float aperture = 0;
    float near = 0.01f;
    float far = 10000;
};

//
// Application state
//
struct yscene {
    // scene data
    yobj::scene* oscn = nullptr;
    ygltf::scene_group* gscn = nullptr;

    // view camera
    ycamera* view_cam = nullptr;

    // view/scene selection
    yobj::camera* ocam = nullptr;
    ygltf::node* gcam = nullptr;

    // time and animation
    float time = 0;
    ym::vec2f time_range = {0, 0};
    bool animate = false;

    // filenames
    std::string filename;
    std::string imfilename;
    std::string outfilename;

    // load options
    bool split_shapes = false;

    // render
    int resolution = 0;
    float exposure = 0, gamma = 2.2f;
    ym::tonemap_type tonemap = ym::tonemap_type::srgb;
    ym::vec4f background = {0, 0, 0, 0};

    // lighting
    bool camera_lights = false;
    ym::vec3f amb = {0, 0, 0};

    // ui
    bool interactive = false;
    bool scene_updated = false;

    // shade
    bool wireframe = false, edges = false;
    bool alpha_cutout = true;
    yshade_state* shstate = nullptr;

    // simulation
    ysym::simulation_params simulation_params;
    int simulation_nframes = 1000;
    ysym::scene* simulation_scene = nullptr;

    // interactive simulation
    std::vector<ym::frame3f> simulation_initial_state;

    // trace
    ytrace::render_params trace_params;
    bool trace_save_progressive = false;
    int trace_block_size = 32;
    int trace_batch_size = 16;
    int trace_nthreads = 0;
    ytrace::scene* trace_scene = nullptr;

    // interactive trace
    int trace_cur_sample = 0, trace_cur_block = 0;
    uint trace_texture_id = 0;
    int trace_preview_width = 0, trace_preview_height = 0;
    int trace_blocks_per_update = 8;
    ym::image<ym::vec4f> trace_hdr, trace_preview;
    std::vector<ym::vec4i> trace_blocks;

    // editing support
    void* selection = nullptr;

    ~yscene() {
        if (oscn) delete oscn;
        if (gscn) delete gscn;
        if (view_cam) delete view_cam;
        if (shstate) delete shstate;
        if (simulation_scene) ysym::free_scene(simulation_scene);
        if (trace_scene) ytrace::free_scene(trace_scene);
    }
};

// ---------------------------------------------------------------------------
// SCENE LOADING
// ---------------------------------------------------------------------------

// loads a scene either OBJ or GLTF
inline bool load_scene(
    yscene* scn, const std::string& filename, bool instances, bool radius) {
    // load scene
    auto ext = yu::path::get_extension(filename);
    if (ext == ".obj") {
        auto err = std::string();
        scn->oscn = yobj::load_scene(filename, true, true, true, true, &err);
        if (!scn->oscn) {
            yu::logging::log_msg(
                yu::logging::log_level::error, "yscene", "cannot load scene");
            yu::logging::log_msg(yu::logging::log_level::error, "ygltf", err);
            return false;
        }
        yobj::add_normals(scn->oscn);
        yobj::add_texture_data(scn->oscn);
        if (radius) yobj::add_radius(scn->oscn, 0.001f);
        yobj::add_tangent_space(scn->oscn);
        if (instances) { yobj::add_instances(scn->oscn); }
        yobj::add_names(scn->oscn);
        if (!scn->oscn->cameras.empty()) {
            auto cam = scn->oscn->cameras[0];
            scn->view_cam = new ycamera();
            scn->view_cam->frame = ym::to_frame(ym::mat4f(cam->xform()));
            scn->view_cam->yfov = cam->yfov;
            scn->view_cam->aspect = cam->aspect;
            scn->view_cam->focus = cam->focus;
            scn->view_cam->aperture = cam->aperture;
        }
        if (scn->split_shapes) yobj::split_shapes(scn->oscn);
    } else if (ext == ".gltf" || ext == ".glb") {
        auto err = std::string();
        scn->gscn = ygltf::load_scenes(filename, true, true, &err);
        if (!scn->gscn) {
            yu::logging::log_msg(
                yu::logging::log_level::error, "yscene", "cannot load scene");
            yu::logging::log_msg(yu::logging::log_level::error, "ygltf", err);
            return false;
        }
        ygltf::add_normals(scn->gscn);
        ygltf::add_texture_data(scn->gscn);
        if (radius) ygltf::add_radius(scn->gscn, 0.001f);
        ygltf::add_tangent_space(scn->gscn);
        if (instances) {
            ygltf::add_nodes(scn->gscn);
            ygltf::add_scene(scn->gscn);
        }
        ygltf::add_names(scn->gscn);
        scn->time_range = ygltf::get_animation_bounds(scn->gscn);
        if (!scn->gscn->default_scene && !scn->gscn->scenes.empty()) {
            scn->gscn->default_scene = scn->gscn->scenes[0];
        }
        if (!ygltf::get_camera_nodes(scn->gscn->default_scene).empty()) {
            auto cam = ygltf::get_camera_nodes(scn->gscn->default_scene)[0];
            scn->view_cam = new ycamera();
            scn->view_cam->frame = ym::to_frame(ym::mat4f(cam->xform()));
            scn->view_cam->yfov = cam->cam->yfov;
            scn->view_cam->aspect = cam->cam->aspect;
            scn->view_cam->near = cam->cam->near;
            scn->view_cam->far = cam->cam->far;
            scn->view_cam->focus = cam->cam->focus;
            scn->view_cam->aperture = cam->cam->aperture;
        }
        if (scn->split_shapes) ygltf::split_shapes(scn->gscn);
    }

#if 1
    if (!scn->view_cam) {
        auto bbox = (scn->oscn) ?
                        ym::bbox3f{yobj::compute_scene_bounds(scn->oscn)} :
                        ym::bbox3f{ygltf::compute_scene_bounds(scn->gscn)};
        auto center = ym::center(bbox);
        auto bbox_size = ym::diagonal(bbox);
        auto bbox_msize = length(bbox_size) / 2;
        // auto bbox_msize =
        //     ym::max(bbox_size[0], ym::max(bbox_size[1], bbox_size[2]));
        scn->view_cam = new ycamera();
        auto camera_dir = ym::vec3f{1.5f, 0.8f, 1.5f};
        auto from = camera_dir * bbox_msize + center;
        auto to = center;
        auto up = ym::vec3f{0, 1, 0};
        scn->view_cam->frame = ym::lookat_frame3(from, to, up);
        scn->view_cam->aspect = 16.0f / 9.0f;
        scn->view_cam->yfov = 2 * atanf(0.5f);
        scn->view_cam->aperture = 0;
        scn->view_cam->focus = ym::length(to - from);
    }
#else
    if (!scn->view_cam) {
        auto bbox = (scn->oscn) ?
                        ym::bbox3f{yobj::compute_scene_bounds(scn->oscn)} :
                        ym::bbox3f{ygltf::compute_scene_bounds(scn->gscn)};
        auto center = ym::center(bbox);
        auto bbox_size = ym::diagonal(bbox);
        scn->view_cam = new ycamera();
        auto from =
            center + ym::vec3f{0, 0, 1} *
                         (bbox_size.z + ym::max(bbox_size.x, bbox_size.y));
        auto to = center;
        auto up = ym::vec3f{0, 1, 0};
        scn->view_cam->frame = ym::lookat_frame3(from, to, up);
        scn->view_cam->aspect = 16.0f / 9.0f;
        scn->view_cam->yfov = 2 * atanf(0.5f);
        scn->view_cam->aperture = 0;
        scn->view_cam->focus = ym::length(to - from);
    }
#endif
    return true;
}

// ---------------------------------------------------------------------------
// DRAWING
// ---------------------------------------------------------------------------

#ifndef YOCTO_NO_OPENGL

//
// Init shading
//
inline void update_shade_lights(yshade_state* st, const yobj::scene* sc) {
    st->lights_pos.clear();
    st->lights_ke.clear();
    st->lights_ltype.clear();

    if (!sc->instances.empty()) {
        for (auto ist : sc->instances) {
            for (auto shp : ist->msh->shapes) {
                if (!shp->mat) continue;
                auto mat = shp->mat;
                if (mat->ke == ym::zero3f) continue;
                for (auto p : shp->points) {
                    if (st->lights_pos.size() >= 16) break;
                    st->lights_pos.push_back(ym::transform_point(
                        ym::mat4f(ist->xform()), shp->pos[p]));
                    st->lights_ke.push_back(mat->ke);
                    st->lights_ltype.push_back(yglu::ltype::point);
                }
            }
        }
    } else {
        for (auto ist : sc->meshes) {
            for (auto shp : ist->shapes) {
                if (!shp->mat) continue;
                auto mat = shp->mat;
                if (mat->ke == ym::zero3f) continue;
                for (auto p : shp->points) {
                    if (st->lights_pos.size() >= 16) break;
                    st->lights_pos.push_back(shp->pos[p]);
                    st->lights_ke.push_back(mat->ke);
                    st->lights_ltype.push_back(yglu::ltype::point);
                }
            }
        }
    }
}

//
// Init shading
//
inline void update_shade_state(const yobj::scene* sc, yshade_state* st) {
    if (!st->prog) yglu::stdshader::make_program(&st->prog, &st->vao);
    st->txt[nullptr] = 0;
    for (auto txt : sc->textures) {
        if (st->txt.find(txt) != st->txt.end()) continue;
        if (txt->hdr) {
            st->txt[txt] = yglu::make_texture(txt->hdr, true, true, true);
        } else if (txt->ldr) {
            st->txt[txt] = yglu::make_texture(txt->ldr, true, true, true);
        } else
            assert(false);
    }
    for (auto mesh : sc->meshes) {
        for (auto shape : mesh->shapes) {
            if (st->vert.find(shape) != st->vert.end()) continue;
            st->vert[shape] = {0, 0, 0, 0};
            st->elem[shape] = {0, 0, 0};
            if (!shape->pos.empty())
                st->vert[shape].pos = yglu::make_buffer((int)shape->pos.size(),
                    3 * sizeof(float), shape->pos.data(), false, false);
            if (!shape->norm.empty())
                st->vert[shape].norm =
                    yglu::make_buffer((int)shape->norm.size(),
                        3 * sizeof(float), shape->norm.data(), false, false);
            if (!shape->texcoord.empty())
                st->vert[shape].texcoord = yglu::make_buffer(
                    (int)shape->texcoord.size(), 2 * sizeof(float),
                    shape->texcoord.data(), false, false);
            if (!shape->color.empty())
                st->vert[shape].color =
                    yglu::make_buffer((int)shape->color.size(),
                        4 * sizeof(float), shape->color.data(), false, false);
            if (!shape->points.empty())
                st->elem[shape].points =
                    yglu::make_buffer((int)shape->points.size(), sizeof(int),
                        shape->points.data(), true, false);
            if (!shape->lines.empty())
                st->elem[shape].lines =
                    yglu::make_buffer((int)shape->lines.size(), 2 * sizeof(int),
                        shape->lines.data(), true, false);
            if (!shape->triangles.empty())
                st->elem[shape].triangles =
                    yglu::make_buffer((int)shape->triangles.size(),
                        3 * sizeof(int), shape->triangles.data(), true, false);
        }
    }
}

//
// Draw a mesh
//
inline void shade_mesh(const yobj::mesh* msh, const yshade_state* st,
    const ym::mat4f& xform, bool edges, bool wireframe, bool cutout) {
    static auto default_material = yobj::material();
    default_material.kd = {0.2f, 0.2f, 0.2f};

    for (auto shp : msh->shapes) {
        yglu::stdshader::begin_shape(st->prog, ym::mat4f(xform));

        auto mat = (shp->mat) ? shp->mat : &default_material;
        yglu::stdshader::set_material_generic(st->prog, mat->ke, mat->kd,
            mat->ks, mat->rs, mat->opacity, st->txt.at(mat->ke_txt),
            st->txt.at(mat->kd_txt), st->txt.at(mat->ks_txt),
            st->txt.at(mat->rs_txt), 0, 0, false, true, cutout);

        auto vert_vbo = st->vert.at(shp);
        yglu::stdshader::set_vert(st->prog, vert_vbo.pos, vert_vbo.norm,
            vert_vbo.texcoord, vert_vbo.color, 0);

        auto velem_vbo = st->elem.at(shp);
        yglu::stdshader::draw_points(
            st->prog, (int)shp->points.size(), velem_vbo.points);
        yglu::stdshader::draw_lines(
            st->prog, (int)shp->lines.size(), velem_vbo.lines);
        yglu::stdshader::draw_triangles(
            st->prog, (int)shp->triangles.size(), velem_vbo.triangles);

        if (edges && !wireframe) {
            assert(yglu::check_error());
            yglu::stdshader::set_material_emission_only(
                st->prog, {0, 0, 0}, mat->opacity, 0, true, cutout);

            assert(yglu::check_error());
            yglu::line_width(2);
            yglu::enable_edges(true);
            yglu::stdshader::draw_points(
                st->prog, (int)shp->points.size(), velem_vbo.points);
            yglu::stdshader::draw_lines(
                st->prog, (int)shp->lines.size(), velem_vbo.lines);
            yglu::stdshader::draw_triangles(
                st->prog, (int)shp->triangles.size(), velem_vbo.triangles);
            yglu::enable_edges(false);
            yglu::line_width(1);
            assert(yglu::check_error());
        }

        yglu::stdshader::end_shape();
    }
}

//
// Display a scene
//
inline void shade_scene(const yobj::scene* sc, yshade_state* st,
    const ycamera* ycam, const yobj::camera* ocam, const ym::vec4f& background,
    float exposure, yglu::tonemap_type tmtype, float gamma, bool wireframe,
    bool edges, bool cutout, bool camera_lights, const ym::vec3f& amb) {
    // update state
    update_shade_state(sc, st);

    // begin frame
    yglu::enable_depth_test(true);
    yglu::enable_culling(false);
    yglu::clear_buffers();

    yglu::enable_wireframe(wireframe);

    ym::mat4f camera_xform, camera_view, camera_proj;
    if (ocam) {
        camera_xform = ym::mat4f(ocam->xform());
        camera_view = ym::inverse(ym::mat4f(ocam->xform()));
        if (ocam->ortho) {
            camera_proj = ym::ortho_mat4(
                ocam->yfov * ocam->aspect, ocam->yfov, 0.01f, 100000.0f);
        } else {
            camera_proj = ym::perspective_mat4(
                ocam->yfov, ocam->aspect, 0.01f, 100000.0f);
        }
    } else {
        camera_xform = ym::to_mat(ycam->frame);
        camera_view = ym::to_mat(ym::inverse(ycam->frame));
        camera_proj = ym::perspective_mat4(
            ycam->yfov, ycam->aspect, ycam->near, ycam->far);
    }

    yglu::stdshader::begin_frame(st->prog, st->vao, camera_lights, exposure,
        tmtype, gamma, camera_xform, camera_view, camera_proj);

    if (!camera_lights) {
        update_shade_lights(st, sc);
        yglu::stdshader::set_lights(st->prog, amb, (int)st->lights_pos.size(),
            st->lights_pos.data(), st->lights_ke.data(),
            st->lights_ltype.data());
    }

    if (!sc->instances.empty()) {
        for (auto ist : sc->instances) {
            shade_mesh(ist->msh, st, ist->xform(), edges, wireframe, cutout);
        }
    } else {
        for (auto msh : sc->meshes) {
            shade_mesh(msh, st, ym::identity_mat4f, edges, wireframe, cutout);
        }
    }

    yglu::stdshader::end_frame();

    yglu::enable_wireframe(false);
}

//
// Init shading
//
inline void update_shade_lights(
    yshade_state* st, const ygltf::scene_group* scn) {
    st->lights_pos.clear();
    st->lights_ke.clear();
    st->lights_ltype.clear();

    auto instances = ygltf::get_mesh_nodes(scn->default_scene);

    if (!instances.empty()) {
        for (auto ist : instances) {
            for (auto shp : ist->msh->shapes) {
                if (!shp->mat) continue;
                auto mat = shp->mat;
                if (mat->emission == ym::zero3f) continue;
                for (auto p : shp->points) {
                    if (st->lights_pos.size() >= 16) break;
                    st->lights_pos.push_back(ym::transform_point(
                        ym::mat4f(ist->xform()), shp->pos[p]));
                    st->lights_ke.push_back(mat->emission);
                    st->lights_ltype.push_back(yglu::ltype::point);
                }
            }
        }
    } else {
        for (auto ist : scn->meshes) {
            for (auto shp : ist->shapes) {
                if (!shp->mat) continue;
                auto mat = shp->mat;
                if (mat->emission == ym::zero3f) continue;
                for (auto p : shp->points) {
                    if (st->lights_pos.size() >= 16) break;
                    st->lights_pos.push_back(shp->pos[p]);
                    st->lights_ke.push_back(mat->emission);
                    st->lights_ltype.push_back(yglu::ltype::point);
                }
            }
        }
    }
}

//
// Init shading
//
inline void update_shade_state(const ygltf::scene_group* sc, yshade_state* st) {
    if (!st->prog) yglu::stdshader::make_program(&st->prog, &st->vao);
    st->txt[nullptr] = 0;
    for (auto txt : sc->textures) {
        if (st->txt.find(txt) != st->txt.end()) continue;
        if (txt->hdr) {
            st->txt[txt] = yglu::make_texture(txt->hdr, true, true, true);
        } else if (txt->ldr) {
            st->txt[txt] = yglu::make_texture(txt->ldr, true, true, true);
        } else
            assert(false);
    }
    for (auto mesh : sc->meshes) {
        for (auto shape : mesh->shapes) {
            if (st->vert.find(shape) != st->vert.end()) continue;
            st->vert[shape] = {0, 0, 0, 0, 0, 0, 0};
            st->elem[shape] = {0, 0, 0};
            if (!shape->pos.empty())
                st->vert[shape].pos = yglu::make_buffer((int)shape->pos.size(),
                    3 * sizeof(float), shape->pos.data(), false, false);
            if (!shape->norm.empty())
                st->vert[shape].norm =
                    yglu::make_buffer((int)shape->norm.size(),
                        3 * sizeof(float), shape->norm.data(), false, false);
            if (!shape->texcoord.empty())
                st->vert[shape].texcoord = yglu::make_buffer(
                    (int)shape->texcoord.size(), 2 * sizeof(float),
                    shape->texcoord.data(), false, false);
            if (!shape->texcoord1.empty())
                st->vert[shape].texcoord1 = yglu::make_buffer(
                    (int)shape->texcoord1.size(), 2 * sizeof(float),
                    shape->texcoord1.data(), false, false);
            if (!shape->color.empty())
                st->vert[shape].color =
                    yglu::make_buffer((int)shape->color.size(),
                        4 * sizeof(float), shape->color.data(), false, false);
            if (!shape->tangsp.empty())
                st->vert[shape].tangsp =
                    yglu::make_buffer((int)shape->tangsp.size(),
                        4 * sizeof(float), shape->tangsp.data(), false, false);
            if (!shape->skin_weights.empty())
                st->vert[shape].skin_weights = yglu::make_buffer(
                    (int)shape->skin_weights.size(), 4 * sizeof(float),
                    shape->skin_weights.data(), false, false);
            if (!shape->skin_joints.empty())
                st->vert[shape].skin_joints = yglu::make_buffer(
                    (int)shape->skin_joints.size(), 4 * sizeof(int),
                    shape->skin_joints.data(), false, false);
            if (!shape->points.empty())
                st->elem[shape].points =
                    yglu::make_buffer((int)shape->points.size(), sizeof(int),
                        shape->points.data(), true, false);
            if (!shape->lines.empty())
                st->elem[shape].lines =
                    yglu::make_buffer((int)shape->lines.size(), 2 * sizeof(int),
                        shape->lines.data(), true, false);
            if (!shape->triangles.empty())
                st->elem[shape].triangles =
                    yglu::make_buffer((int)shape->triangles.size(),
                        3 * sizeof(int), shape->triangles.data(), true, false);
        }
    }
}

//
// Draw a mesh
//
inline void shade_mesh(const ygltf::mesh* msh, const ygltf::skin* sk,
    const std::vector<float>& morph_weights, const yshade_state* st,
    const ym::mat4f& xform, bool edges, bool wireframe, bool cutout) {
    static auto default_material = ygltf::material();
    default_material.metallic_roughness =
        new ygltf::material_metallic_rooughness();
    default_material.metallic_roughness->base = {0.2f, 0.2f, 0.2f};

    auto txt_info =
        [st](ygltf::texture* gtxt,
            const ygltf::texture_info* ginfo) -> yglu::texture_info {
        auto info = yglu::texture_info();
        info.txt_id = 0;
        if (!gtxt) return info;
        info.txt_id = st->txt.at(gtxt);
        if (!ginfo) return info;
        info.wrap_s = (yglu::texture_wrap)ginfo->wrap_s;
        info.wrap_t = (yglu::texture_wrap)ginfo->wrap_t;
        info.filter_min = (yglu::texture_filter)ginfo->filter_min;
        info.filter_mag = (yglu::texture_filter)ginfo->filter_mag;
        return info;
    };

    for (auto shp : msh->shapes) {
        yglu::stdshader::begin_shape(st->prog, xform);

        auto mat = (shp->mat) ? shp->mat : &default_material;
        float op = 1;
        if (mat->specular_glossiness) {
            auto sg = mat->specular_glossiness;
            op = sg->opacity;
            yglu::stdshader::set_material_gltf_specular_glossiness(st->prog,
                mat->emission, sg->diffuse, sg->specular, sg->glossiness,
                sg->opacity,
                txt_info(mat->emission_txt, mat->emission_txt_info),
                txt_info(sg->diffuse_txt, sg->diffuse_txt_info),
                txt_info(sg->specular_txt, sg->specular_txt_info),
                txt_info(mat->normal_txt, mat->normal_txt_info),
                txt_info(mat->occlusion_txt, mat->occlusion_txt_info), false,
                mat->double_sided, cutout);
        } else if (mat->metallic_roughness) {
            auto mr = mat->metallic_roughness;
            op = mr->opacity;
            yglu::stdshader::set_material_gltf_metallic_roughness(st->prog,
                mat->emission, mr->base, mr->metallic, mr->roughness,
                mr->opacity,
                txt_info(mat->emission_txt, mat->emission_txt_info),
                txt_info(mr->base_txt, mr->base_txt_info),
                txt_info(mr->metallic_txt, mr->metallic_txt_info),
                txt_info(mat->normal_txt, mat->normal_txt_info),
                txt_info(mat->occlusion_txt, mat->occlusion_txt_info), false,
                mat->double_sided, cutout);
        } else {
            yglu::stdshader::set_material_emission_only(st->prog, mat->emission,
                1, txt_info(mat->emission_txt, mat->emission_txt_info),
                mat->double_sided, cutout);
        }

        auto vert_vbo = st->vert.at(shp);
        yglu::stdshader::set_vert(st->prog, vert_vbo.pos, vert_vbo.norm,
            vert_vbo.texcoord, vert_vbo.color, vert_vbo.tangsp);
        if (sk) {
            auto skin_xforms = ygltf::get_skin_transforms(sk, xform);
#if 0
            yglu::stdshader::set_vert_gltf_skinning(st->prog,
                vert_vbo.skin_weights, vert_vbo.skin_joints, skin_xforms.size(),
                skin_xforms.data());
#else
            std::vector<ym::vec3f> skinned_pos, skinned_norm;
            // ym::compute_skinning(shp->pos, shp->norm, shp->skin_weights,
            //     shp->skin_joints, skin_xforms, skinned_pos, skinned_norm);
            ym::compute_matrix_skinning(shp->pos, shp->norm, shp->skin_weights,
                shp->skin_joints, skin_xforms, skinned_pos, skinned_norm);
            yglu::update_buffer(vert_vbo.pos, 3 * sizeof(float),
                (int)skinned_pos.size(), skinned_pos.data(), false, true);
            yglu::update_buffer(vert_vbo.norm, 3 * sizeof(float),
                (int)skinned_norm.size(), skinned_norm.data(), false, true);
#endif
        } else if (!morph_weights.empty() && !shp->morph_targets.empty()) {
            std::vector<ym::vec3f> morph_pos, morph_norm;
            std::vector<ym::vec4f> morph_tang;
            ygltf::compute_morphing_deformation(
                shp, morph_weights, morph_pos, morph_norm, morph_tang);
            yglu::update_buffer(vert_vbo.pos, 3 * sizeof(float),
                (int)morph_pos.size(), morph_pos.data(), false, true);
            yglu::update_buffer(vert_vbo.norm, 3 * sizeof(float),
                (int)morph_norm.size(), morph_norm.data(), false, true);
        } else {
            yglu::stdshader::set_vert_skinning_off(st->prog);
        }

        yglu::enable_culling(!mat->double_sided);

        auto velem_vbo = st->elem.at(shp);
        yglu::stdshader::draw_points(
            st->prog, (int)shp->points.size(), velem_vbo.points);
        yglu::stdshader::draw_lines(
            st->prog, (int)shp->lines.size(), velem_vbo.lines);
        yglu::stdshader::draw_triangles(
            st->prog, (int)shp->triangles.size(), velem_vbo.triangles);

        yglu::enable_culling(false);

        if (edges && !wireframe) {
            assert(yglu::check_error());
            yglu::stdshader::set_material_emission_only(
                st->prog, {0, 0, 0}, op, 0, true, cutout);

            assert(yglu::check_error());
            yglu::line_width(2);
            yglu::enable_edges(true);
            yglu::stdshader::draw_points(
                st->prog, (int)shp->points.size(), velem_vbo.points);
            yglu::stdshader::draw_lines(
                st->prog, (int)shp->lines.size(), velem_vbo.lines);
            yglu::stdshader::draw_triangles(
                st->prog, (int)shp->triangles.size(), velem_vbo.triangles);
            yglu::enable_edges(false);
            yglu::line_width(1);
            assert(yglu::check_error());
        }

        yglu::stdshader::end_shape();
    }
}

//
// Display a scene
//
inline void shade_scene(const ygltf::scene_group* scns, yshade_state* st,
    const ycamera* ycam, const ygltf::node* gcam, const ym::vec4f& background,
    float exposure, yglu::tonemap_type tmtype, float gamma, bool wireframe,
    bool edges, bool cutout, bool camera_lights, const ym::vec3f& amb) {
    // update state
    update_shade_state(scns, st);

    // begin frame
    yglu::enable_depth_test(true);
    yglu::enable_culling(false);
    yglu::clear_buffers();

    yglu::enable_wireframe(wireframe);

    ym::mat4f camera_xform, camera_view, camera_proj;
    if (gcam) {
        camera_xform = ym::mat4f(gcam->xform());
        camera_view = ym::inverse(ym::mat4f(gcam->xform()));
        if (gcam->cam->ortho) {
            auto near = (gcam->cam->near) ? gcam->cam->near : 0.001f;
            auto far = (gcam->cam->far) ? gcam->cam->far : 10000;
            camera_proj = ym::ortho_mat4(gcam->cam->yfov * gcam->cam->aspect,
                gcam->cam->yfov, near, far);
        } else {
            auto near = (gcam->cam->near) ? gcam->cam->near : 0.001f;
            if (gcam->cam->far) {
                camera_proj = ym::perspective_mat4(
                    gcam->cam->yfov, gcam->cam->aspect, near, gcam->cam->far);
            } else {
                camera_proj = ym::perspective_mat4(
                    gcam->cam->yfov, gcam->cam->aspect, 0.01f);
            }
        }
    } else {
        camera_xform = ym::to_mat(ycam->frame);
        camera_view = ym::to_mat(ym::inverse(ycam->frame));
        auto near = (ycam->near) ? ycam->near : 0.001f;
        if (ycam->far) {
            camera_proj =
                ym::perspective_mat4(ycam->yfov, ycam->aspect, near, ycam->far);
        } else {
            camera_proj = ym::perspective_mat4(ycam->yfov, ycam->aspect, near);
        }
    }

    yglu::stdshader::begin_frame(st->prog, st->vao, camera_lights, exposure,
        tmtype, gamma, camera_xform, camera_view, camera_proj);

    if (!camera_lights) {
        update_shade_lights(st, scns);
        yglu::stdshader::set_lights(st->prog, amb, (int)st->lights_pos.size(),
            st->lights_pos.data(), st->lights_ke.data(),
            st->lights_ltype.data());
    }

    auto instances = ygltf::get_mesh_nodes(scns->default_scene);

    if (!instances.empty()) {
        for (auto ist : instances) {
            shade_mesh(ist->msh, ist->skn, ist->morph_weights, st, ist->xform(),
                edges, wireframe, cutout);
        }
    } else {
        for (auto msh : scns->meshes) {
            shade_mesh(msh, nullptr, {}, st, ym::identity_mat4f, edges,
                wireframe, cutout);
        }
    }

    yglu::stdshader::end_frame();

    yglu::enable_wireframe(false);
}

#endif

// ---------------------------------------------------------------------------
// INTERACTIVE FUNCTIONS
// ---------------------------------------------------------------------------

void save_screenshot(ygui::window* win, const std::string& imfilename) {
    if (yu::path::get_extension(imfilename) != ".png") {
        printf("supports only png screenshots");
        return;
    }

    auto wh = ym::vec2i{0, 0};
    auto pixels = ygui::get_screenshot(win, wh);
    yimg::save_image(
        imfilename, wh[0], wh[1], 4, (unsigned char*)pixels.data());
}

void draw_scene(ygui::window* win) {
    auto scn = (yscene*)ygui::get_user_pointer(win);
    auto window_size = ygui::get_window_size(win);
    auto aspect = (float)window_size[0] / (float)window_size[1];
    scn->view_cam->aspect = aspect;
    if (scn->ocam) scn->ocam->aspect = aspect;
    if (scn->gcam) scn->gcam->cam->aspect = aspect;
    if (scn->oscn) {
        shade_scene(scn->oscn, scn->shstate, scn->view_cam, scn->ocam,
            scn->background, scn->exposure, (yglu::tonemap_type)scn->tonemap,
            scn->gamma, scn->wireframe, scn->edges, scn->alpha_cutout,
            scn->camera_lights, scn->amb);
    } else {
        shade_scene(scn->gscn, scn->shstate, scn->view_cam, scn->gcam,
            scn->background, scn->exposure, (yglu::tonemap_type)scn->tonemap,
            scn->gamma, scn->wireframe, scn->edges, scn->alpha_cutout,
            scn->camera_lights, scn->amb);
    }
}

void draw_tree_widgets(ygui::window* win, const std::string& lbl,
    yobj::camera* cam, void** selection) {
    ygui::tree_leaf_widget(win, lbl + cam->name, selection, cam);
}

void draw_tree_widgets(ygui::window* win, const std::string& lbl,
    yobj::texture* txt, void** selection) {
    ygui::tree_leaf_widget(win, lbl + txt->path, selection, txt);
}

void draw_tree_widgets(ygui::window* win, const std::string& lbl,
    yobj::material* mat, void** selection) {
    if (ygui::tree_begin_widget(win, lbl + mat->name, selection, mat)) {
        if (mat->ke_txt) draw_tree_widgets(win, "ke: ", mat->ke_txt, selection);
        if (mat->kd_txt) draw_tree_widgets(win, "ke: ", mat->kd_txt, selection);
        if (mat->ks_txt) draw_tree_widgets(win, "ke: ", mat->ks_txt, selection);
        ygui::tree_end_widget(win);
    }
}

void draw_tree_widgets(ygui::window* win, const std::string& lbl,
    yobj::shape* shp, void** selection) {
    if (ygui::tree_begin_widget(win, lbl + shp->name, selection, shp)) {
        if (shp->mat) draw_tree_widgets(win, "mat: ", shp->mat, selection);
        ygui::tree_end_widget(win);
    }
}

void draw_tree_widgets(ygui::window* win, const std::string& lbl,
    yobj::mesh* msh, void** selection) {
    if (ygui::tree_begin_widget(win, lbl + msh->name, selection, msh)) {
        auto sid = 0;
        for (auto shp : msh->shapes) {
            draw_tree_widgets(
                win, "sh" + std::to_string(sid++) + ": ", shp, selection);
        }
        ygui::tree_end_widget(win);
    }
}

void draw_tree_widgets(ygui::window* win, const std::string& lbl,
    yobj::instance* ist, void** selection) {
    if (ygui::tree_begin_widget(win, lbl + ist->name, selection, ist)) {
        if (ist->msh) draw_tree_widgets(win, "mesh: ", ist->msh, selection);
        ygui::tree_end_widget(win);
    }
}

void draw_tree_widgets(ygui::window* win, const std::string& lbl,
    yobj::scene* oscn, void** selection) {
    if (ygui::tree_begin_widget(win, lbl + "cameras")) {
        for (auto cam : oscn->cameras)
            draw_tree_widgets(win, "", cam, selection);
        ygui::tree_end_widget(win);
    }

    if (ygui::tree_begin_widget(win, lbl + "meshes")) {
        for (auto msh : oscn->meshes)
            draw_tree_widgets(win, "", msh, selection);
        ygui::tree_end_widget(win);
    }

    if (ygui::tree_begin_widget(win, lbl + "instances")) {
        for (auto ist : oscn->instances)
            draw_tree_widgets(win, "", ist, selection);
        ygui::tree_end_widget(win);
    }

    if (ygui::tree_begin_widget(win, lbl + "materials")) {
        for (auto mat : oscn->materials)
            draw_tree_widgets(win, "", mat, selection);
        ygui::tree_end_widget(win);
    }

    if (ygui::tree_begin_widget(win, lbl + "textures")) {
        for (auto txt : oscn->textures)
            draw_tree_widgets(win, "", txt, selection);
        ygui::tree_end_widget(win);
    }
}

void draw_elem_widgets(ygui::window* win, yobj::scene* oscn, yobj::texture* txt,
    void** selection, const yshade_state* state) {
    if (selection && *selection != txt) return;

    ygui::separator_widget(win);
    ygui::label_widget(win, "path", txt->path);
    auto size = yu::string::format("%d x %d @ 4 %s", txt->width(),
        txt->height(), (txt->ldr) ? "byte" : "float");
    ygui::label_widget(win, "size", size);
    ygui::image_widget(
        win, state->txt.at(txt), {128, 128}, {txt->width(), txt->height()});
}

void draw_elem_widgets(ygui::window* win, yobj::scene* oscn,
    yobj::material* mat, void** selection, const yshade_state* state) {
    if (selection && *selection != mat) return;

    auto txt_names = std::vector<std::pair<std::string, yobj::texture*>>{
        {"<none>", nullptr}};
    for (auto txt : oscn->textures) txt_names.push_back({txt->path, txt});

    ygui::separator_widget(win);
    ygui::label_widget(win, "name", mat->name);
    ygui::slider_widget(win, "ke", &mat->ke, 0, 1000);
    ygui::slider_widget(win, "kd", &mat->kd, 0, 1);
    ygui::slider_widget(win, "ks", &mat->ks, 0, 1);
    ygui::slider_widget(win, "kt", &mat->kt, 0, 1);
    ygui::slider_widget(win, "rs", &mat->rs, 0, 1);
    ygui::combo_widget(win, "ke_txt", &mat->ke_txt, txt_names);
    ygui::combo_widget(win, "kd_txt", &mat->kd_txt, txt_names);
    ygui::combo_widget(win, "ks_txt", &mat->ks_txt, txt_names);
    ygui::combo_widget(win, "kt_txt", &mat->kt_txt, txt_names);
}

void draw_elem_widgets(ygui::window* win, yobj::scene* oscn, yobj::shape* shp,
    void** selection, const yshade_state* state) {
    if (selection && *selection != shp) return;

    auto mat_names = std::vector<std::pair<std::string, yobj::material*>>{
        {"<none>", nullptr}};
    for (auto mat : oscn->materials) mat_names.push_back({mat->name, mat});

    ygui::separator_widget(win);
    ygui::label_widget(win, "name", shp->name);
    ygui::combo_widget(win, "material", &shp->mat, mat_names);
    ygui::label_widget(win, "verts", (int)shp->pos.size());
    if (!shp->triangles.empty())
        ygui::label_widget(win, "triangles", (int)shp->triangles.size());
    if (!shp->lines.empty())
        ygui::label_widget(win, "lines", (int)shp->lines.size());
    if (!shp->points.empty())
        ygui::label_widget(win, "points", (int)shp->points.size());
    if (!shp->tetras.empty())
        ygui::label_widget(win, "tetras", (int)shp->tetras.size());
}

void draw_elem_widgets(ygui::window* win, yobj::scene* oscn, yobj::mesh* msh,
    void** selection, const yshade_state* state) {
    if (selection && *selection != msh) return;
    ygui::separator_widget(win);
    ygui::label_widget(win, "name", msh->name);
    for (auto shp : msh->shapes) {
        draw_elem_widgets(win, oscn, shp, nullptr, state);
    }
}

void draw_elem_widgets(ygui::window* win, yobj::scene* oscn, yobj::camera* cam,
    void** selection, const yshade_state* state) {
    if (selection && *selection != cam) return;
    ygui::separator_widget(win);
    ygui::label_widget(win, "name", cam->name);
    ygui::slider_widget(win, "translation", &cam->translation, -10, 10);
    ygui::slider_widget(win, "rotation", &cam->rotation);
    ygui::slider_widget(win, "matrix", &cam->matrix, -10, 10);
    ygui::slider_widget(win, "yfov", &cam->yfov, 0.1, 4);
    ygui::slider_widget(win, "aspect", &cam->aspect, 0.1, 4);
    ygui::slider_widget(win, "focus", &cam->focus, 0.01, 10);
    ygui::slider_widget(win, "aperture", &cam->aperture, 0, 1);
}

void draw_elem_widgets(ygui::window* win, yobj::scene* oscn,
    yobj::instance* ist, void** selection, const yshade_state* state) {
    if (selection && *selection != ist) return;

    auto msh_names =
        std::vector<std::pair<std::string, yobj::mesh*>>{{"<none>", nullptr}};
    for (auto msh : oscn->meshes) msh_names.push_back({msh->name, msh});

    ygui::separator_widget(win);
    ygui::label_widget(win, "name", ist->name);
    ygui::slider_widget(win, "translation", &ist->translation, -10, 10);
    ygui::slider_widget(win, "rotation", &ist->rotation);
    ygui::slider_widget(win, "scale", &ist->scale, 0.01, 10);
    ygui::slider_widget(win, "matrix", &ist->matrix, -10, 10);
    ygui::combo_widget(win, "mesh", &ist->msh, msh_names);
}

void draw_elem_widgets(ygui::window* win, yobj::scene* oscn, void** selection,
    const yshade_state* state) {
    for (auto cam : oscn->cameras) {
        draw_elem_widgets(win, oscn, cam, selection, state);
    }

    for (auto msh : oscn->meshes) {
        draw_elem_widgets(win, oscn, msh, selection, state);
        for (auto shp : msh->shapes) {
            draw_elem_widgets(win, oscn, shp, selection, state);
        }
    }

    for (auto ist : oscn->instances) {
        draw_elem_widgets(win, oscn, ist, selection, state);
    }

    for (auto mat : oscn->materials) {
        draw_elem_widgets(win, oscn, mat, selection, state);
    }

    for (auto txt : oscn->textures) {
        draw_elem_widgets(win, oscn, txt, selection, state);
    }
}

void draw_scene_widgets(ygui::window* win, yobj::scene* oscn, void** selection,
    const yshade_state* state) {
    ygui::scroll_region_begin_widget(win, "model", 240, false);
    draw_tree_widgets(win, "", oscn, selection);
    ygui::scroll_region_end_widget(win);
    draw_elem_widgets(win, oscn, selection, state);
}

void draw_tree_widgets(ygui::window* win, const std::string& lbl,
    ygltf::camera* cam, void** selection) {
    ygui::tree_leaf_widget(win, lbl + cam->name, selection, cam);
}

void draw_tree_widgets(ygui::window* win, const std::string& lbl,
    ygltf::texture* txt, void** selection) {
    ygui::tree_leaf_widget(win, lbl + txt->path, selection, txt);
}

void draw_tree_widgets(ygui::window* win, const std::string& lbl,
    ygltf::material_metallic_rooughness* mat, void** selection) {
    if (mat->base_txt) draw_tree_widgets(win, "kb: ", mat->base_txt, selection);
    if (mat->metallic_txt)
        draw_tree_widgets(win, "km: ", mat->metallic_txt, selection);
}

void draw_tree_widgets(ygui::window* win, const std::string& lbl,
    ygltf::material_specular_glossiness* mat, void** selection) {
    if (mat->diffuse_txt)
        draw_tree_widgets(win, "kd: ", mat->diffuse_txt, selection);
    if (mat->specular_txt)
        draw_tree_widgets(win, "ks: ", mat->specular_txt, selection);
}

void draw_tree_widgets(ygui::window* win, const std::string& lbl,
    ygltf::material* mat, void** selection) {
    if (ygui::tree_begin_widget(win, lbl + mat->name, selection, mat)) {
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
        ygui::tree_end_widget(win);
    }
}

void draw_tree_widgets(ygui::window* win, const std::string& lbl,
    ygltf::shape* shp, void** selection) {
    if (ygui::tree_begin_widget(win, lbl + shp->name, selection, shp)) {
        if (shp->mat) draw_tree_widgets(win, "mat: ", shp->mat, selection);
        ygui::tree_end_widget(win);
    }
}

void draw_tree_widgets(ygui::window* win, const std::string& lbl,
    ygltf::mesh* msh, void** selection) {
    if (ygui::tree_begin_widget(win, lbl + msh->name, selection, msh)) {
        auto sid = 0;
        for (auto shp : msh->shapes) {
            draw_tree_widgets(
                win, "sh" + std::to_string(sid++) + ": ", shp, selection);
        }
        ygui::tree_end_widget(win);
    }
}

void draw_tree_widgets(ygui::window* win, const std::string& lbl,
    ygltf::node* node, void** selection) {
    if (ygui::tree_begin_widget(win, lbl + node->name, selection, node)) {
        if (node->msh) draw_tree_widgets(win, "mesh: ", node->msh, selection);
        if (node->cam) draw_tree_widgets(win, "cam: ", node->cam, selection);
        for (auto child : node->children)
            draw_tree_widgets(win, "", child, selection);
        ygui::tree_end_widget(win);
    }
}

void draw_tree_widgets(ygui::window* win, const std::string& lbl,
    ygltf::scene* scn, void** selection) {
    if (ygui::tree_begin_widget(win, lbl + scn->name, selection, scn)) {
        for (auto node : scn->nodes)
            draw_tree_widgets(win, "", node, selection);
        ygui::tree_end_widget(win);
    }
}

void draw_tree_widgets(ygui::window* win, const std::string& lbl,
    ygltf::animation* anim, void** selection) {
    if (ygui::tree_begin_widget(win, lbl, selection, anim)) {
        auto sid = 0;
        for (auto node : anim->nodes) {
            draw_tree_widgets(
                win, "node " + std::to_string(sid++) + ": ", node, selection);
        }
        ygui::tree_end_widget(win);
    }
}

void draw_tree_widgets(ygui::window* win, const std::string& lbl,
    ygltf::animation_group* anims, void** selection) {
    if (ygui::tree_begin_widget(win, lbl + anims->name, selection, anims)) {
        auto sid = 0;
        for (auto anim : anims->animations) {
            draw_tree_widgets(
                win, "anim " + std::to_string(sid++), anim, selection);
        }
        ygui::tree_end_widget(win);
    }
}

void draw_tree_widgets(ygui::window* win, const std::string& lbl,
    ygltf::skin* skin, void** selection) {
    if (ygui::tree_begin_widget(win, lbl + skin->name, selection, skin)) {
        draw_tree_widgets(win, "skeleton", skin->root, selection);
        ygui::tree_end_widget(win);
    }
}

void draw_tree_widgets(ygui::window* win, const std::string& lbl,
    ygltf::scene_group* oscn, void** selection) {
    if (ygui::tree_begin_widget(win, lbl + "scenes")) {
        for (auto scn : oscn->scenes)
            draw_tree_widgets(win, "", scn, selection);
        ygui::tree_end_widget(win);
    }

    if (ygui::tree_begin_widget(win, lbl + "cameras")) {
        for (auto cam : oscn->cameras)
            draw_tree_widgets(win, "", cam, selection);
        ygui::tree_end_widget(win);
    }

    if (ygui::tree_begin_widget(win, lbl + "meshes")) {
        for (auto msh : oscn->meshes)
            draw_tree_widgets(win, "", msh, selection);
        ygui::tree_end_widget(win);
    }

    if (ygui::tree_begin_widget(win, lbl + "nodes")) {
        for (auto node : oscn->nodes) {
            if (!node->parent) draw_tree_widgets(win, "", node, selection);
        }
        ygui::tree_end_widget(win);
    }

    if (ygui::tree_begin_widget(win, lbl + "materials")) {
        for (auto mat : oscn->materials)
            draw_tree_widgets(win, "", mat, selection);
        ygui::tree_end_widget(win);
    }

    if (ygui::tree_begin_widget(win, lbl + "textures")) {
        for (auto txt : oscn->textures)
            draw_tree_widgets(win, "", txt, selection);
        ygui::tree_end_widget(win);
    }

    if (ygui::tree_begin_widget(win, lbl + "animations")) {
        for (auto anim : oscn->animations)
            draw_tree_widgets(win, "", anim, selection);
        ygui::tree_end_widget(win);
    }

    if (ygui::tree_begin_widget(win, lbl + "skins")) {
        for (auto skin : oscn->skins)
            draw_tree_widgets(win, "", skin, selection);
        ygui::tree_end_widget(win);
    }
}

void draw_elem_widgets(ygui::window* win, ygltf::scene_group* gscn,
    ygltf::texture* txt, void** selection, const yshade_state* state) {
    if (selection && *selection != txt) return;

    ygui::separator_widget(win);
    ygui::label_widget(win, "path", txt->path);
    auto str = yu::string::format("%d x %d @ 4 %s", txt->width(), txt->height(),
        (txt->ldr) ? "byte" : "float");
    ygui::label_widget(win, "size", str);
    if (state->txt.find(txt) != state->txt.end()) {
        ygui::image_widget(
            win, state->txt.at(txt), {128, 128}, {txt->width(), txt->height()});
    }
}

void draw_elem_widgets(ygui::window* win, ygltf::scene_group* gscn,
    ygltf::texture_info* info, void** selection, const yshade_state* state) {
    static const auto wrap_names =
        std::vector<std::pair<std::string, ygltf::texture_wrap>>{
            {"repeat", ygltf::texture_wrap::repeat},
            {"clamp", ygltf::texture_wrap::clamp},
            {"mirror", ygltf::texture_wrap::mirror}};
    static const auto filter_mag_names =
        std::vector<std::pair<std::string, ygltf::texture_filter>>{
            {"linear", ygltf::texture_filter::linear},
            {"nearest", ygltf::texture_filter::nearest},
        };
    static const auto filter_min_names =
        std::vector<std::pair<std::string, ygltf::texture_filter>>{
            {"linear", ygltf::texture_filter::linear},
            {"nearest", ygltf::texture_filter::nearest},
            {"linear_mipmap_linear",
                ygltf::texture_filter::linear_mipmap_linear},
            {"linear_mipmap_nearest",
                ygltf::texture_filter::linear_mipmap_nearest},
            {"nearest_mipmap_linear",
                ygltf::texture_filter::nearest_mipmap_linear},
            {"nearest_mipmap_nearest",
                ygltf::texture_filter::nearest_mipmap_nearest},
        };

    if (selection && *selection != info) return;

    ygui::indent_begin_widgets(win);

    ygui::combo_widget(win, "wrap s", &info->wrap_s, wrap_names);
    ygui::combo_widget(win, "wrap t", &info->wrap_t, wrap_names);
    ygui::combo_widget(win, "filter mag", &info->filter_mag, filter_mag_names);
    ygui::combo_widget(win, "filter min", &info->filter_min, filter_min_names);

    ygui::indent_end_widgets(win);
}

void draw_elem_widgets(ygui::window* win, ygltf::scene_group* gscn,
    ygltf::material_metallic_rooughness* mat, void** selection,
    const yshade_state* state) {
    auto txt_names = std::vector<std::pair<std::string, ygltf::texture*>>{
        {"<none>", nullptr}};
    for (auto txt : gscn->textures) txt_names.push_back({txt->path, txt});

    ygui::slider_widget(win, "mr base", &mat->base, 0, 1);
    ygui::slider_widget(win, "mr opacity", &mat->opacity, 0, 1);
    ygui::slider_widget(win, "mr roughness", &mat->roughness, 0, 1);
    ygui::combo_widget(win, "mr base txt", &mat->base_txt, txt_names);
    if (mat->base_txt && mat->base_txt_info)
        draw_elem_widgets(win, gscn, mat->base_txt_info, nullptr, state);
    ygui::combo_widget(win, "mr metallic txt", &mat->metallic_txt, txt_names);
    if (mat->metallic_txt && mat->metallic_txt_info)
        draw_elem_widgets(win, gscn, mat->metallic_txt_info, nullptr, state);
}

void draw_elem_widgets(ygui::window* win, ygltf::scene_group* gscn,
    ygltf::material_specular_glossiness* mat, void** selection,
    const yshade_state* state) {
    auto txt_names = std::vector<std::pair<std::string, ygltf::texture*>>{
        {"<none>", nullptr}};
    for (auto txt : gscn->textures) txt_names.push_back({txt->path, txt});

    ygui::slider_widget(win, "sg diffuse", &mat->diffuse, 0, 1);
    ygui::slider_widget(win, "sg opacity", &mat->opacity, 0, 1);
    ygui::slider_widget(win, "sg specular", &mat->specular, 0, 1);
    ygui::slider_widget(win, "sg glossiness", &mat->glossiness, 0, 1);
    ygui::combo_widget(win, "sg diffuse txt", &mat->diffuse_txt, txt_names);
    if (mat->diffuse_txt && mat->diffuse_txt_info)
        draw_elem_widgets(win, gscn, mat->diffuse_txt_info, nullptr, state);
    ygui::combo_widget(win, "sg specular txt", &mat->specular_txt, txt_names);
    if (mat->specular_txt && mat->specular_txt_info)
        draw_elem_widgets(win, gscn, mat->specular_txt_info, nullptr, state);
}

void draw_elem_widgets(ygui::window* win, ygltf::scene_group* gscn,
    ygltf::material* mat, void** selection, const yshade_state* state) {
    if (selection && *selection != mat) return;

    auto txt_names = std::vector<std::pair<std::string, ygltf::texture*>>{
        {"<none>", nullptr}};
    for (auto txt : gscn->textures) txt_names.push_back({txt->path, txt});

    ygui::separator_widget(win);
    ygui::label_widget(win, "name", mat->name);
    ygui::slider_widget(win, "emission", &mat->emission, 0, 1000);
    ygui::combo_widget(win, "emission txt", &mat->emission_txt, txt_names);
    if (mat->emission_txt && mat->emission_txt_info)
        draw_elem_widgets(win, gscn, mat->emission_txt_info, nullptr, state);
    ygui::combo_widget(win, "normal txt", &mat->normal_txt, txt_names);
    if (mat->normal_txt && mat->normal_txt_info)
        draw_elem_widgets(win, gscn, mat->normal_txt_info, nullptr, state);
    ygui::combo_widget(win, "occlusion txt", &mat->occlusion_txt, txt_names);
    if (mat->occlusion_txt && mat->occlusion_txt_info)
        draw_elem_widgets(win, gscn, mat->occlusion_txt_info, nullptr, state);
    if (mat->metallic_roughness) {
        ygui::separator_widget(win);
        draw_elem_widgets(win, gscn, mat->metallic_roughness, nullptr, state);
    } else {
        if (ygui::button_widget(win, "add metallic roughness")) {
            mat->metallic_roughness = new ygltf::material_metallic_rooughness();
        }
    }
    if (mat->specular_glossiness) {
        ygui::separator_widget(win);
        draw_elem_widgets(win, gscn, mat->specular_glossiness, nullptr, state);
    } else {
        if (ygui::button_widget(win, "add specular glossiness")) {
            mat->specular_glossiness =
                new ygltf::material_specular_glossiness();
        }
    }
}

void draw_elem_widgets(ygui::window* win, ygltf::scene_group* gscn,
    ygltf::shape* shp, void** selection, const yshade_state* state) {
    if (selection && *selection != shp) return;

    auto mat_names = std::vector<std::pair<std::string, ygltf::material*>>{
        {"<none>", nullptr}};
    for (auto mat : gscn->materials) mat_names.push_back({mat->name, mat});

    ygui::separator_widget(win);
    ygui::label_widget(win, "name", shp->name);
    ygui::combo_widget(win, "material", &shp->mat, mat_names);
    ygui::label_widget(win, "verts", (int)shp->pos.size());
    if (!shp->triangles.empty())
        ygui::label_widget(win, "triangles", (int)shp->triangles.size());
    if (!shp->lines.empty())
        ygui::label_widget(win, "lines", (int)shp->lines.size());
    if (!shp->points.empty())
        ygui::label_widget(win, "points", (int)shp->points.size());
}

void draw_elem_widgets(ygui::window* win, ygltf::scene_group* gscn,
    ygltf::mesh* msh, void** selection, const yshade_state* state) {
    if (selection && *selection != msh) return;
    ygui::separator_widget(win);
    ygui::label_widget(win, "name", msh->name);
    for (auto shp : msh->shapes) {
        draw_elem_widgets(win, gscn, shp, nullptr, state);
    }
}

void draw_elem_widgets(ygui::window* win, ygltf::scene_group* gscn,
    ygltf::camera* cam, void** selection, const yshade_state* state) {
    if (selection && *selection != cam) return;
    ygui::separator_widget(win);
    ygui::label_widget(win, "name", cam->name);
    ygui::slider_widget(win, "yfov", &cam->yfov, 0.1, 4);
    ygui::slider_widget(win, "aspect", &cam->aspect, 0.1, 4);
    ygui::slider_widget(win, "focus", &cam->focus, 0.01, 10);
    ygui::slider_widget(win, "aperture", &cam->aperture, 0, 1);
}

void draw_elem_widgets(ygui::window* win, ygltf::scene_group* gscn,
    ygltf::animation* anim, void** selection, const yshade_state* state) {
    if (selection && *selection != anim) return;
    auto interp_name = std::map<ygltf::animation_interpolation, std::string>{
        {ygltf::animation_interpolation::step, "step"},
        {ygltf::animation_interpolation::linear, "linear"},
        {ygltf::animation_interpolation::catmull_rom, "catmull-rom spline"},
        {ygltf::animation_interpolation::cubic, "cubic spline"},
    };
    auto node_names =
        std::vector<std::pair<std::string, ygltf::node*>>{{"<none>", nullptr}};
    for (auto node : gscn->nodes) node_names.push_back({node->name, node});
    ygui::separator_widget(win);
    ygui::label_widget(win, "interpolation", interp_name.at(anim->interp));
    ygui::label_widget(win, "num times", (int)anim->time.size());
    ygui::label_widget(win, "num translation", (int)anim->translation.size());
    ygui::label_widget(win, "num rotation", (int)anim->rotation.size());
    ygui::label_widget(win, "num scale", (int)anim->scale.size());
    auto selected_node = (ygltf::node*)nullptr;
    ygui::list_widget(win, "", &selected_node, node_names);
}

void draw_elem_widgets(ygui::window* win, ygltf::scene_group* gscn,
    ygltf::animation_group* anims, void** selection,
    const yshade_state* state) {
    if (selection && *selection != anims) return;
    ygui::separator_widget(win);
    ygui::label_widget(win, "name", anims->name);
    for (auto anim : anims->animations) {
        draw_elem_widgets(win, gscn, anim, nullptr, state);
    }
}

void draw_elem_widgets(ygui::window* win, ygltf::scene_group* gscn,
    ygltf::skin* skin, void** selection, const yshade_state* state) {
    if (selection && *selection != skin) return;
    ygui::separator_widget(win);
    ygui::label_widget(win, "name", skin->name);
}

void draw_elem_widgets(ygui::window* win, ygltf::scene_group* gscn,
    ygltf::node* node, void** selection, const yshade_state* state) {
    if (selection && *selection != node) return;

    auto msh_names =
        std::vector<std::pair<std::string, ygltf::mesh*>>{{"<none>", nullptr}};
    for (auto msh : gscn->meshes) msh_names.push_back({msh->name, msh});
    auto cam_names = std::vector<std::pair<std::string, ygltf::camera*>>{
        {"<none>", nullptr}};
    for (auto cam : gscn->cameras) cam_names.push_back({cam->name, cam});

    ygui::separator_widget(win);
    ygui::label_widget(win, "name", node->name);
    ygui::combo_widget(win, "mesh", &node->msh, msh_names);
    ygui::combo_widget(win, "cam", &node->cam, cam_names);
    ygui::slider_widget(win, "translation", &node->translation, -10, 10);
    ygui::slider_widget(win, "rotation", &node->rotation);
    ygui::slider_widget(win, "scale", &node->scale, 0.01, 10);
    ygui::slider_widget(win, "matrix", &node->matrix, -10, 10);
}

void draw_elem_widgets(ygui::window* win, ygltf::scene_group* gscn,
    void** selection, const yshade_state* state) {
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

    for (auto anim : gscn->animations) {
        draw_elem_widgets(win, gscn, anim, selection, state);
    }

    for (auto skin : gscn->skins) {
        draw_elem_widgets(win, gscn, skin, selection, state);
    }
}

template <typename T>
inline void remove(std::vector<T*>& v, T* vv) {
    auto pos = std::find(v.begin(), v.end(), vv);
    if (pos != v.end()) v.erase(pos);
}

inline void draw_edit_widgets(ygui::window* win, ygltf::scene_group* gscn,
    void** selection, const yshade_state* state) {
    static auto shape_names =
        std::vector<std::pair<std::string, int>>{{"cube", 0}, {"sphere", 1}};
    static auto shape_type = 0;
    static char txt_filename[1024] = "grid.png";

    auto selected_node = (ygltf::node*)nullptr;
    auto selected_scene = (ygltf::scene*)nullptr;
    auto selected_mesh = (ygltf::mesh*)nullptr;
    auto selected_cam = (ygltf::camera*)nullptr;
    auto selected_txt = (ygltf::texture*)nullptr;
    auto selected_mat = (ygltf::material*)nullptr;

    if (*selection) {
        for (auto ptr : gscn->nodes)
            if (ptr == *selection) selected_node = ptr;
        for (auto ptr : gscn->scenes)
            if (ptr == *selection) selected_scene = ptr;
        for (auto ptr : gscn->meshes)
            if (ptr == *selection) selected_mesh = ptr;
        for (auto ptr : gscn->cameras)
            if (ptr == *selection) selected_cam = ptr;
        for (auto ptr : gscn->materials)
            if (ptr == *selection) selected_mat = ptr;
        for (auto ptr : gscn->textures)
            if (ptr == *selection) selected_txt = ptr;
    }

    static auto auto_parent = true;
    ygui::checkbox_widget(win, "set parent from selection", &auto_parent);

    if (ygui::button_widget(win, "add mesh")) {
        static auto count = 0;
        auto mesh = new ygltf::mesh();
        mesh->name = "<new mesh " + std::to_string(count++) + ">";
        auto shp = new ygltf::shape();
        mesh->shapes.push_back(shp);
        shp->name = "<new shape " + std::to_string(count - 1) + ">";
        switch (shape_type) {
            case 0: {
                ym::make_uvcube(
                    1, 1, shp->triangles, shp->pos, shp->norm, shp->texcoord);
            } break;
            case 1: {
                ym::make_uvsphere(
                    1, 1, shp->triangles, shp->pos, shp->norm, shp->texcoord);
            } break;
        }
        if (auto_parent && selected_node) selected_node->msh = mesh;
        gscn->meshes.push_back(mesh);
        *selection = mesh;
    }
    ygui::combo_widget(win, "shape type", &shape_type, shape_names);

    if (ygui::button_widget(win, "add camera")) {
        static auto count = 0;
        auto cam = new ygltf::camera();
        cam->name = "<new camera " + std::to_string(count++) + ">";
        if (auto_parent && selected_node) selected_node->cam = cam;
        gscn->cameras.push_back(cam);
        *selection = cam;
    }

    if (ygui::button_widget(win, "add node")) {
        static auto count = 0;
        auto node = new ygltf::node();
        node->name = "<new node " + std::to_string(count++) + ">";
        if (auto_parent && selected_node)
            selected_node->children.push_back(node);
        gscn->nodes.push_back(node);
        *selection = node;
    }

    if (ygui::button_widget(win, "add texture")) {
        static auto count = 0;
        auto txt = new ygltf::texture();
        txt->name = "<new texture " + std::to_string(count++) + ">";
        txt->path = txt_filename;
        auto scn = (yscene*)ygui::get_user_pointer(win);
        auto dirname = yu::path::get_dirname(scn->filename);
        try {
            auto ext = yu::path::get_extension(dirname + txt->path);
            if (ext == ".hdr") {
                txt->hdr = yimg::load_image4f(dirname + txt->path);
            } else {
                txt->ldr = yimg::load_image4b(dirname + txt->path);
            }
        } catch (...) {
            txt->ldr = new ym::image4b(1, 1, {255, 255, 255, 255});
        }
        gscn->textures.push_back(txt);
        *selection = txt;
    }
    ygui::text_widget(win, "texture", txt_filename, sizeof(txt_filename));

    if (ygui::button_widget(win, "delete")) {
        if (selected_cam) {
            for (auto node : gscn->nodes)
                if (node->cam == selected_cam) node->cam = nullptr;
            remove(gscn->cameras, selected_cam);
            delete selected_cam;
            *selection = nullptr;
        }
        if (selected_mesh) {
            for (auto node : gscn->nodes)
                if (node->msh == selected_mesh) node->cam = nullptr;
            remove(gscn->meshes, selected_mesh);
            delete selected_mesh;
            *selection = nullptr;
        }
        if (selected_mat) {
            for (auto mesh : gscn->meshes)
                for (auto shp : mesh->shapes)
                    if (shp->mat == selected_mat) shp->mat = nullptr;
            remove(gscn->materials, selected_mat);
            delete selected_mat;
            *selection = nullptr;
        }
        if (selected_txt) {
            for (auto mat : gscn->materials) {
                if (mat->emission_txt == selected_txt)
                    mat->emission_txt = nullptr;
                if (mat->normal_txt == selected_txt) mat->normal_txt = nullptr;
                if (mat->occlusion_txt == selected_txt)
                    mat->occlusion_txt = nullptr;
                if (mat->metallic_roughness) {
                    if (mat->metallic_roughness->base_txt == selected_txt)
                        mat->metallic_roughness->base_txt = nullptr;
                    if (mat->metallic_roughness->metallic_txt == selected_txt)
                        mat->metallic_roughness->metallic_txt = nullptr;
                }
                if (mat->specular_glossiness) {
                    if (mat->specular_glossiness->diffuse_txt == selected_txt)
                        mat->specular_glossiness->diffuse_txt = nullptr;
                    if (mat->specular_glossiness->specular_txt == selected_txt)
                        mat->specular_glossiness->specular_txt = nullptr;
                }
            }
            remove(gscn->textures, selected_txt);
            delete selected_txt;
            *selection = nullptr;
        }
    }
}

void draw_scene_widgets(ygui::window* win, ygltf::scene_group* gscn,
    void** selection, const yshade_state* state) {
    ygui::scroll_region_begin_widget(win, "model", 240, false);
    draw_tree_widgets(win, "", gscn, selection);
    ygui::scroll_region_end_widget(win);
    draw_elem_widgets(win, gscn, selection, state);
}

#ifndef YSCENE_SIMULATE_HACK
void simulate_step(yscene* scene) {}
#else
void simulate_step(yscene* scene);
#endif

void draw_widgets(ygui::window* win) {
    static auto tmtype_names = std::vector<std::pair<std::string, int>>{
        {"none", (int)ym::tonemap_type::none},
        {"srgb", (int)ym::tonemap_type::srgb},
        {"gamma", (int)ym::tonemap_type::gamma},
        {"filmic", (int)ym::tonemap_type::filmic}};
    static auto rtype_names = std::vector<std::pair<std::string, int>>{
        {"uniform", (int)ytrace::rng_type::uniform},
        {"stratified", (int)ytrace::rng_type::stratified},
        {"cmjs", (int)ytrace::rng_type::cmjs}};
    static auto stype_names = std::vector<std::pair<std::string, int>>{
        {"eye", (int)ytrace::shader_type::eyelight},
        {"direct", (int)ytrace::shader_type::direct},
        {"direct_ao", (int)ytrace::shader_type::direct_ao},
        {"path", (int)ytrace::shader_type::pathtrace}};

    auto scn = (yscene*)ygui::get_user_pointer(win);
    if (ygui::begin_widgets(win, "yview")) {
        ygui::label_widget(win, "scene", scn->filename);
        if (scn->time_range != ym::zero2f && !scn->simulation_scene &&
            !scn->trace_scene) {
            ygui::slider_widget(
                win, "time", &scn->time, scn->time_range.x, scn->time_range.y);
            ygui::checkbox_widget(win, "play animation", &scn->animate);
        }
        if (scn->simulation_scene) {
            ygui::slider_widget(
                win, "time", &scn->time, scn->time_range.x, scn->time_range.y);
            if (ygui::button_widget(win, "start")) scn->animate = true;
            ygui::continue_line_widgets(win);
            if (ygui::button_widget(win, "stop")) scn->animate = false;
            ygui::continue_line_widgets(win);
            if (ygui::button_widget(win, "step")) { simulate_step(scn); }
        }
        if (scn->trace_scene) {
            ygui::label_widget(win, "s", scn->trace_cur_sample);
            ygui::slider_widget(
                win, "samples", &scn->trace_params.nsamples, 0, 1000000, 1);
            if (ygui::combo_widget(win, "shader type",
                    (int*)&scn->trace_params.stype, stype_names))
                scn->scene_updated = true;
            if (ygui::combo_widget(win, "random type",
                    (int*)&scn->trace_params.rtype, rtype_names))
                scn->scene_updated = true;
        }
        if (scn->oscn) {
            auto camera_names =
                std::vector<std::pair<std::string, yobj::camera*>>{
                    {"<view>", nullptr}};
            for (auto cam : scn->oscn->cameras)
                camera_names.push_back({cam->name, cam});
            ygui::combo_widget(win, "camera", &scn->ocam, camera_names);
        } else if (scn->gscn && scn->gscn->default_scene) {
            auto camera_names =
                std::vector<std::pair<std::string, ygltf::node*>>{
                    {"<view>", nullptr}};
            for (auto cam : ygltf::get_camera_nodes(scn->gscn->default_scene)) {
                camera_names.push_back({cam->name, cam});
            }
            ygui::combo_widget(win, "camera", &scn->gcam, camera_names);
        }
        if (!scn->trace_scene) {
            ygui::checkbox_widget(win, "eyelight", &scn->camera_lights);
            ygui::continue_line_widgets(win);
            ygui::checkbox_widget(win, "wire", &scn->wireframe);
            ygui::continue_line_widgets(win);
            ygui::checkbox_widget(win, "edges", &scn->edges);
            ygui::continue_line_widgets(win);
            ygui::checkbox_widget(win, "cutout", &scn->alpha_cutout);
        }
        ygui::slider_widget(win, "hdr exposure", &scn->exposure, -20, 20, 1);
        ygui::slider_widget(win, "hdr gamma", &scn->gamma, 0.1, 5, 0.1);
        ygui::combo_widget(
            win, "hdr tonemap", (int*)&scn->tonemap, tmtype_names);
        if (ygui::collapsing_header_widget(win, "view cam")) {
            auto cam = scn->view_cam;
            ygui::slider_widget(win, "yfov", &cam->yfov, 0.1, 5);
            ygui::slider_widget(win, "aspect", &cam->aspect, 0.1, 3);
            ygui::slider_widget(win, "focus", &cam->focus, 0.1, 1000);
            ygui::slider_widget(win, "near", &cam->near, 0.01, 1);
            ygui::slider_widget(win, "far", &cam->far, 1, 10000);
        }
        if (ygui::collapsing_header_widget(win, "edit")) {
            if (scn->gscn)
                draw_edit_widgets(
                    win, scn->gscn, &scn->selection, scn->shstate);
        }
        if (ygui::collapsing_header_widget(win, "inspect")) {
            if (scn->oscn)
                draw_scene_widgets(
                    win, scn->oscn, &scn->selection, scn->shstate);
            if (scn->gscn)
                draw_scene_widgets(
                    win, scn->gscn, &scn->selection, scn->shstate);
        }
    }
    ygui::end_widgets(win);
}

// ---------------------------------------------------------------------------
// COMMAND LINE
// ---------------------------------------------------------------------------

void parse_cmdline(yscene* scene, int argc, char** argv, const char* help,
    bool simulation, bool trace) {
    static auto rtype_names = std::vector<std::pair<std::string, int>>{
        {"uniform", (int)ytrace::rng_type::uniform},
        {"stratified", (int)ytrace::rng_type::stratified},
        {"cmjs", (int)ytrace::rng_type::cmjs}};
    static auto stype_names = std::vector<std::pair<std::string, int>>{
        {"path", (int)ytrace::shader_type::pathtrace},
        {"eye", (int)ytrace::shader_type::eyelight},
        {"direct", (int)ytrace::shader_type::direct},
        {"direct_ao", (int)ytrace::shader_type::direct_ao}};
    static auto tmtype_names = std::vector<std::pair<std::string, int>>{
        {"none", (int)ym::tonemap_type::none},
        {"srgb", (int)ym::tonemap_type::srgb},
        {"gamma", (int)ym::tonemap_type::gamma},
        {"filmic", (int)ym::tonemap_type::filmic}};

    // parse command line
    auto parser = yu::cmdline::make_parser(argc, argv, help);

#ifndef YOCTO_NO_OPENGL
    if (simulation || trace) {
        scene->interactive =
            parse_flag(parser, "--interactive", "-i", "interactive mode");
    }
#endif

    // simulation
    if (simulation) {
        scene->simulation_params.dt =
            parse_optf(parser, "--delta_time", "-dt", "delta time", 1 / 60.0f);
        scene->simulation_nframes =
            parse_opti(parser, "--nframes", "-n", "number of frames", 1000);
        scene->outfilename = parse_opts(
            parser, "--output-scene", "-O", "output filename", "out.%04d.obj");
    }

    // trace
    if (trace) {
        scene->trace_params.camera_id = 0;
        scene->trace_save_progressive = parse_flag(
            parser, "--save-progressive", "", "save progressive images");

        scene->trace_params.rtype =
            (ytrace::rng_type)parse_opte(parser, "--random", "", "random type",
                (int)ytrace::rng_type::stratified, rtype_names);
        scene->trace_params.stype = (ytrace::shader_type)parse_opte(parser,
            "--shader", "-S", "path estimator type",
            (int)ytrace::shader_type::pathtrace, stype_names);
        scene->trace_params.envmap_invisible =
            parse_flag(parser, "--envmap-invisible", "", "envmap invisible");
        scene->trace_nthreads = parse_opti(
            parser, "--threads", "-t", "number of threads [0 for default]", 0);
        scene->trace_block_size =
            parse_opti(parser, "--block-size", "", "block size", 32);
        scene->trace_batch_size =
            parse_opti(parser, "--batch-size", "", "batch size", 16);
        scene->trace_params.nsamples =
            parse_opti(parser, "--samples", "-s", "image samples", 256);
    }

    // render
    scene->exposure =
        parse_optf(parser, "--exposure", "-e", "hdr image exposure", 0);
    scene->gamma = parse_optf(parser, "--gamma", "-g", "hdr image gamma", 2.2f);
    scene->tonemap = (ym::tonemap_type)parse_opte(parser, "--tonemap", "-t",
        "hdr tonemap output", (int)ym::tonemap_type::srgb, tmtype_names);
    scene->resolution =
        parse_opti(parser, "--resolution", "-r", "image resolution", 720);

    auto amb = parse_optf(parser, "--ambient", "", "ambient factor", 0);

    scene->camera_lights = parse_flag(
        parser, "--camera-lights", "-c", "enable camera lights", false);

    scene->amb = {amb, amb, amb};
    scene->trace_params.amb = {amb, amb, amb};
    if (scene->camera_lights) {
        scene->trace_params.stype = ytrace::shader_type::eyelight;
    }

    // load options
    scene->split_shapes = parse_flag(
        parser, "--split-shapes", "", "split meshes into single shapes", false);

    // params
    scene->imfilename =
        parse_opts(parser, "--output-image", "-o", "image filename", "out.hdr");
    scene->filename = parse_args(parser, "scene", "scene filename", "", true);

    // check parsing
    check_parser(parser);
}

//
// Logging
//
inline void set_default_loggers() {
    auto loggers = yu::logging::get_default_loggers();
    loggers->push_back(yu::logging::make_stdout_logger());
    loggers->push_back(
        make_file_logger("yocto.log", true, yu::logging::log_level::verbose));
}

// ---------------------------------------------------------------------------
// INTERACTION LOOP
// ---------------------------------------------------------------------------

//
// callbacks
//
using init_fn = void (*)(ygui::window* win);
using draw_fn = void (*)(ygui::window* win);
using update_fn = bool (*)(yscene* scn);

//
// init draw with shading
//
void shade_init(ygui::window* win) {
    auto scn = (yscene*)ygui::get_user_pointer(win);
    scn->shstate = new yshade_state();
    if (scn->oscn) {
        update_shade_state(scn->oscn, scn->shstate);
    } else if (scn->gscn) {
        update_shade_state(scn->gscn, scn->shstate);
    }
}

//
// draw with shading
//
void shade_draw(ygui::window* win) {
    auto scn = (yscene*)ygui::get_user_pointer(win);
    if (scn->gscn) {
        ygltf::update_animated_transforms(scn->gscn, scn->time);
        ygltf::update_transforms(scn->gscn);
    }
    draw_scene(win);
    draw_widgets(win);
    ygui::swap_buffers(win);
    if (scn->shstate->lights_pos.empty()) scn->camera_lights = true;
}

//
// run ui loop
//
void run_ui(yscene* scn, int w, int h, const std::string& title, init_fn init,
    draw_fn draw, update_fn update) {
    // window
    auto win = ygui::init_window(w, h, title, scn);
    ygui::set_callbacks(win, nullptr, nullptr, draw);

    // window values
    int mouse_button = 0;
    ym::vec2f mouse_pos, mouse_last;

    // load textures
    init(win);

    // init widget
    ygui::init_widgets(win);

    // loop
    while (!ygui::should_close(win)) {
        mouse_last = mouse_pos;
        mouse_pos = ygui::get_mouse_posf(win);
        mouse_button = ygui::get_mouse_button(win);

        ygui::set_window_title(win, ("yshade | " + scn->filename));

        // handle mouse
        if (!scn->gcam && !scn->ocam && mouse_button &&
            mouse_pos != mouse_last && !ygui::get_widget_active(win)) {
            auto dolly = 0.0f;
            auto pan = ym::zero2f;
            auto rotate = ym::zero2f;
            switch (mouse_button) {
                case 1: rotate = (mouse_pos - mouse_last) / 100.0f; break;
                case 2: dolly = (mouse_pos[0] - mouse_last[0]) / 100.0f; break;
                case 3: pan = (mouse_pos - mouse_last) / 100.0f; break;
                default: break;
            }

            ym::turntable(
                scn->view_cam->frame, scn->view_cam->focus, rotate, dolly, pan);
            scn->scene_updated = true;
        }

        // draw
        draw(win);

        // update
        auto updated = update(scn);

        // check for screenshot
        //        if (scn->no_ui) {
        //            save_screenshot(win, scn->imfilename);
        //            break;
        //        }

        // event hadling
        if (updated)
            ygui::poll_events(win);
        else
            ygui::wait_events(win);
    }

    ygui::clear_widgets(win);
    ygui::clear_window(win);
}

#endif
