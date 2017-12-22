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

#define YGL_OPENGL 1
#include "../yocto/yocto_gl.h"
using namespace ygl;

// ---------------------------------------------------------------------------
// SCENE (OBJ or GLTF) AND APPLICATION PARAMETERS
// ---------------------------------------------------------------------------

//
// OpenGL shape vbo
//
struct shape_vbo {
    gl_vertex_buffer pos = {};
    gl_vertex_buffer norm = {};
    gl_vertex_buffer texcoord = {};
    gl_vertex_buffer texcoord1 = {};
    gl_vertex_buffer color = {};
    gl_vertex_buffer tangsp = {};
    gl_element_buffer points = {};
    gl_element_buffer lines = {};
    gl_element_buffer triangles = {};
    gl_element_buffer quads = {};
    gl_element_buffer edges = {};
};

//
// OpenGL state
//
struct shade_state {
    // shade state
    gl_stdsurface_program prog = {};
    unordered_map<texture*, gl_texture> txt;
    unordered_map<shape*, shape_vbo> vbo;

    // lights
    vector<vec3f> lights_pos;
    vector<vec3f> lights_ke;
    vector<gl_ltype> lights_ltype;
};

//
// Application state
//
struct app_state {
    // scene data
    scene* scn = nullptr;

    // camera selection
    camera* scam = nullptr;

    // filenames
    string filename;
    string imfilename;
    string outfilename;

    // render
    int resolution = 0;
    float exposure = 0, gamma = 2.2f;
    bool filmic = false;
    vec4f background = {0, 0, 0, 0};

    // lighting
    bool camera_lights = false;
    vec3f amb = {0, 0, 0};

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
        if (scn) delete scn;
    }
};

// Init shading
inline void update_shade_lights(shade_state* st, const scene* scn) {
    st->lights_pos.clear();
    st->lights_ke.clear();
    st->lights_ltype.clear();

    if (!scn->instances.empty()) {
        for (auto ist : scn->instances) {
            auto shp = ist->shp;
            if (!shp->mat) continue;
            if (shp->mat->ke == zero3f) continue;
            if (!shp->points.empty()) {
                for (auto p : shp->points) {
                    if (st->lights_pos.size() >= 16) break;
                    st->lights_pos += transform_point(ist->frame, shp->pos[p]);
                    st->lights_ke += shp->mat->ke;
                    st->lights_ltype += gl_ltype::point;
                }
            } else {
                auto bbox = make_bbox(shp->pos.size(), shp->pos.data());
                auto pos = bbox_center(bbox);
                auto area = 0.0f;
                for (auto l : shp->lines)
                    area += line_length(shp->pos[l.x], shp->pos[l.y]);
                for (auto t : shp->triangles)
                    area += triangle_area(
                        shp->pos[t.x], shp->pos[t.y], shp->pos[t.z]);
                for (auto t : shp->quads)
                    area += quad_area(shp->pos[t.x], shp->pos[t.y],
                        shp->pos[t.z], shp->pos[t.w]);
                auto ke = shp->mat->ke * area;
                if (st->lights_pos.size() < 16) {
                    st->lights_pos += transform_point(ist->frame, pos);
                    st->lights_ke += ke;
                    st->lights_ltype += gl_ltype::point;
                }
            }
        }
    } else {
        for (auto shp : scn->shapes) {
            if (!shp->mat) continue;
            if (shp->mat->ke == zero3f) continue;
            for (auto p : shp->points) {
                if (st->lights_pos.size() >= 16) break;
                st->lights_pos.push_back(shp->pos[p]);
                st->lights_ke.push_back(shp->mat->ke);
                st->lights_ltype.push_back(gl_ltype::point);
            }
        }
    }
}

// Init shading
inline void update_shade_state(const scene* scn, shade_state* st) {
    if (!is_program_valid(st->prog)) st->prog = make_stdsurface_program();
    st->txt[nullptr] = {};
    for (auto txt : scn->textures) {
        if (st->txt.find(txt) != st->txt.end()) continue;
        if (txt->hdr) {
            st->txt[txt] = make_texture(txt->hdr, true, true, true);
        } else if (txt->ldr) {
            st->txt[txt] = make_texture(txt->ldr, true, true, true);
        } else
            assert(false);
    }
    for (auto shp : scn->shapes) {
        if (st->vbo.find(shp) != st->vbo.end()) continue;
        st->vbo[shp] = shape_vbo();
        if (!shp->pos.empty()) st->vbo[shp].pos = make_vertex_buffer(shp->pos);
        if (!shp->norm.empty())
            st->vbo[shp].norm = make_vertex_buffer(shp->norm);
        if (!shp->texcoord.empty())
            st->vbo[shp].texcoord = make_vertex_buffer(shp->texcoord);
        if (!shp->color.empty())
            st->vbo[shp].color = make_vertex_buffer(shp->color);
        if (!shp->tangsp.empty())
            st->vbo[shp].tangsp = make_vertex_buffer(shp->tangsp);
        if (!shp->points.empty())
            st->vbo[shp].points = make_element_buffer(shp->points);
        if (!shp->lines.empty())
            st->vbo[shp].lines = make_element_buffer(shp->lines);
        if (!shp->triangles.empty()) {
            st->vbo[shp].triangles = make_element_buffer(shp->triangles);
        }
        if (!shp->quads.empty()) {
            auto triangles = convert_quads_to_triangles(shp->quads);
            st->vbo[shp].quads = make_element_buffer(triangles);
        }
        if (!shp->triangles.empty() || !shp->quads.empty() ||
            !shp->quads_pos.empty()) {
            auto edges = get_edges(
                shp->lines, shp->triangles, shp->quads + shp->quads_pos);
            st->vbo[shp].edges = make_element_buffer(edges);
        }
    }
}

// Draw a shape
inline void shade_shape(shape* shp, shade_state* st, const mat4f& xform,
    bool highlighted, bool edges, bool wireframe, bool cutout) {
    static auto default_material = material();
    default_material.kd = {0.2f, 0.2f, 0.2f};

    begin_stdsurface_shape(st->prog, mat4f(xform));

    auto etype = gl_etype::triangle;
    if (!shp->lines.empty()) etype = gl_etype::line;
    if (!shp->points.empty()) etype = gl_etype::point;

    set_stdsurface_highlight(
        st->prog, (highlighted) ? vec4f{1, 1, 0, 1} : zero4f);

    auto txt = [&st](texture_info& info) -> gl_texture_info {
        if (!info.txt) return {};
        return st->txt.at(info.txt);
    };

    auto mat = (shp->mat) ? shp->mat : &default_material;
    set_stdsurface_material(st->prog, mat->mtype, etype, mat->ke, mat->kd,
        mat->ks, mat->rs, mat->op, txt(mat->ke_txt), txt(mat->kd_txt),
        txt(mat->ks_txt), txt(mat->rs_txt), txt(mat->norm_txt),
        txt(mat->occ_txt), false, mat->double_sided, cutout);

    auto& vbo = st->vbo.at(shp);
    set_stdsurface_vert(
        st->prog, vbo.pos, vbo.norm, vbo.texcoord, vbo.color, vbo.tangsp);

    draw_elems(vbo.points);
    draw_elems(vbo.lines);
    draw_elems(vbo.triangles);
    draw_elems(vbo.quads);

    if (edges && !wireframe) {
        assert(gl_check_error());
        set_stdsurface_material(st->prog, material_type::specular_roughness,
            etype, zero3f, zero3f, zero3f, 0.5f, mat->op, {}, {}, {}, {}, {},
            {}, true, mat->double_sided, cutout);

        assert(gl_check_error());
        gl_line_width(2);
        gl_enable_edges(true);
        draw_elems(vbo.edges);
        gl_enable_edges(false);
        gl_line_width(1);
        assert(gl_check_error());
    }

    end_stdsurface_shape(st->prog);
}

// Display a scene
inline void shade_scene(const scene* scn, shade_state* st, const camera* cam,
    void* selection, const vec4f& background, float exposure, float gamma,
    bool filmic, bool wireframe, bool edges, bool cutout, bool camera_lights,
    const vec3f& amb) {
    // update state
    update_shade_state(scn, st);

    // begin frame
    gl_enable_depth_test(true);
    gl_enable_culling(false);

    gl_enable_wireframe(wireframe);

    mat4f camera_xform, camera_view, camera_proj;
    camera_xform = to_mat4f(cam->frame);
    camera_view = to_mat4f(inverse(cam->frame));
    camera_proj =
        perspective_mat4f(cam->yfov, cam->aspect, cam->near, cam->far);

    begin_stdsurface_frame(st->prog, camera_lights, exposure, gamma, filmic,
        camera_xform, camera_view, camera_proj);

    if (!camera_lights) {
        update_shade_lights(st, scn);
        set_stdsurface_lights(st->prog, amb, (int)st->lights_pos.size(),
            st->lights_pos.data(), st->lights_ke.data(),
            st->lights_ltype.data());
    }

    if (!scn->instances.empty()) {
        for (auto ist : scn->instances) {
            shade_shape(ist->shp, st, ist->xform(),
                (ist == selection || ist->shp == selection), edges, wireframe,
                cutout);
        }
    } else {
        for (auto shp : scn->shapes) {
            shade_shape(shp, st, identity_mat4f, shp == selection, edges,
                wireframe, cutout);
        }
    }

    end_stdsurface_frame(st->prog);
    gl_enable_wireframe(false);
}

// draw with shading
inline void draw(gl_window* win) {
    auto app = (app_state*)get_user_pointer(win);

    auto window_size = get_window_size(win);
    auto framebuffer_size = get_framebuffer_size(win);
    gl_set_viewport(framebuffer_size);
    auto aspect = (float)window_size.x / (float)window_size.y;
    app->scam->aspect = aspect;

    gl_clear_buffers();
    shade_scene(app->scn, app->shstate, app->scam, app->selection,
        app->background, app->exposure, app->gamma, app->filmic, app->wireframe,
        app->edges, app->alpha_cutout, app->camera_lights, app->amb);

    if (begin_widgets(win, "yview")) {
        draw_label_widget(win, "scene", app->filename);
        draw_camera_widget(win, "camera", app->scn, app->scam);
        draw_value_widget(win, "wire", app->wireframe);
        draw_continue_widget(win);
        draw_value_widget(win, "edges", app->edges);
        draw_continue_widget(win);
        draw_value_widget(win, "cutout", app->alpha_cutout);
        draw_continue_widget(win);
        draw_value_widget(win, "fps", app->navigation_fps);
        draw_tonemap_widgets(win, "", app->exposure, app->gamma, app->filmic);
        draw_scene_widgets(
            win, "scene", app->scn, app->selection, app->shstate->txt);
    }
    end_widgets(win);

    swap_buffers(win);
    if (app->shstate->lights_pos.empty()) app->camera_lights = true;
}

// scene update
bool update(app_state* st) { return false; }

//
// run ui loop
//
inline void run_ui(app_state* app, int w, int h, const string& title) {
    // window
    auto win = make_window(w, h, title, app);
    set_window_callbacks(win, nullptr, nullptr, draw);

    // window values
    int mouse_button = 0;
    vec2f mouse_pos, mouse_last;

    // load textures and vbos
    app->shstate = new shade_state();
    update_shade_state(app->scn, app->shstate);

    // init widget
    init_widgets(win);

    // loop
    while (!should_close(win)) {
        mouse_last = mouse_pos;
        mouse_pos = get_mouse_posf(win);
        mouse_button = get_mouse_button(win);

        set_window_title(win, ("yview | " + app->filename));

        // handle mouse and keyboard for navigation
        if (mouse_button && !get_widget_active(win)) {
            if (app->navigation_fps) {
                auto dolly = 0.0f;
                auto pan = zero2f;
                auto rotate = zero2f;
                switch (mouse_button) {
                    case 1: rotate = (mouse_pos - mouse_last) / 100.0f; break;
                    case 2:
                        dolly = (mouse_pos[0] - mouse_last[0]) / 100.0f;
                        break;
                    case 3: pan = (mouse_pos - mouse_last) / 100.0f; break;
                    default: break;
                }
                camera_fps(app->scam->frame, {0, 0, 0}, rotate);
            } else {
                auto dolly = 0.0f;
                auto pan = zero2f;
                auto rotate = zero2f;
                switch (mouse_button) {
                    case 1: rotate = (mouse_pos - mouse_last) / 100.0f; break;
                    case 2:
                        dolly = (mouse_pos[0] - mouse_last[0]) / 100.0f;
                        break;
                    case 3: pan = (mouse_pos - mouse_last) / 100.0f; break;
                    default: break;
                }

                camera_turntable(
                    app->scam->frame, app->scam->focus, rotate, dolly, pan);
            }
            app->scene_updated = true;
        }

        // handle keytboard for navigation
        if (!get_widget_active(win) && app->navigation_fps) {
            auto transl = zero3f;
            if (get_key(win, 'a')) transl.x -= 1;
            if (get_key(win, 'd')) transl.x += 1;
            if (get_key(win, 's')) transl.z += 1;
            if (get_key(win, 'w')) transl.z -= 1;
            if (get_key(win, 'e')) transl.y += 1;
            if (get_key(win, 'q')) transl.y -= 1;
            if (transl != zero3f) {
                camera_fps(app->scam->frame, transl, {0, 0});
                app->scene_updated = true;
            }
        }

        // draw
        draw(win);

        // update
        update(app);

        // check for screenshot
        //        if (scn->no_ui) {
        //            save_screenshot(win, scn->imfilename);
        //            break;
        //        }

        // event hadling
        poll_events(win);
    }

    clear_window(win);
}

// ---------------------------------------------------------------------------
// MAIN
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    // create empty scene
    auto app = new app_state();

    // parse command line
    auto parser = make_parser(argc, argv, "yview", "views scenes inteactively");
    app->exposure =
        parse_opt(parser, "--exposure", "-e", "hdr image exposure", 0.0f);
    app->gamma = parse_opt(parser, "--gamma", "-g", "hdr image gamma", 2.2f);
    app->filmic = parse_flag(parser, "--filmic", "-F", "hdr filmic output");
    app->resolution =
        parse_opt(parser, "--resolution", "-r", "image resolution", 540);
    auto amb = parse_opt(parser, "--ambient", "", "ambient factor", 0.0f);
    app->amb = {amb, amb, amb};
    app->camera_lights =
        parse_flag(parser, "--camera-lights", "-c", "enable camera lights");
    auto log_filename = parse_opt(parser, "--log", "", "log to disk", ""s);
    if (log_filename != "") add_file_stream(log_filename, true);
    auto preserve_quads =
        parse_flag(parser, "--preserve-quads", "-q", "preserve quads on load");
    auto preserve_facevarying = parse_flag(
        parser, "--preserve-facevarying", "-f", "preserve facevarying on load");
    app->imfilename =
        parse_opt(parser, "--output-image", "-o", "image filename", "out.hdr"s);
    app->filename = parse_arg(parser, "scene", "scene filename", ""s);
    if (should_exit(parser)) {
        printf("%s\n", get_usage(parser).c_str());
        exit(1);
    }

    // scene loading
    log_info("loading scene {}", app->filename);
    try {
        auto opts = load_options();
        opts.preserve_quads = preserve_quads;
        opts.preserve_facevarying = preserve_facevarying;
        app->scn = load_scene(app->filename, opts);
    } catch (exception e) { log_fatal("cannot load scene {}", app->filename); }

    // tesselate input shapes
    tesselate_shapes(app->scn);

    // add missing data
    add_elements(app->scn);
    app->scam = app->scn->cameras[0];

    // light
    update_lights(app->scn, true);

    // run ui
    auto width = (int)round(app->scam->aspect * app->resolution);
    auto height = app->resolution;
    run_ui(app, width, height, "yview");

    // clear
    delete app;

    // done
    return 0;
}
