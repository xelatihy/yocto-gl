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

#include "yapp.h"

#include "../yocto/yocto_glu.h"
#include "../yocto/yocto_img.h"
#include "../yocto/yocto_math.h"

struct state {
    // params
    yapp::params* pars = nullptr;

    // scene
    yapp::scene* scene = nullptr;
    ybvh::scene* scene_bvh = nullptr;
    ysym::scene* rigid_scene = nullptr;

    // animation state
    bool simulating = false;
    int frame = 0;
    std::vector<ym::frame3f> initial_state;
    bool show_debug = false;

    // shade state
    yapp::shade_state* shst = nullptr;

    // widgets
    void* widget_ctx = nullptr;

    // clear
    ~state() {
        if (pars) delete pars;
        if (scene) delete scene;
        if (scene_bvh) ybvh::free_scene(scene_bvh);
        if (rigid_scene) ysym::free_scene(rigid_scene);
        if (shst) yapp::free_shade_state(shst);
    }
};

const int hud_width = 256;

void save_screenshot(yglu::ui::window* win, const std::string& imfilename) {
    if (ycmd::get_extension(imfilename) != ".png") {
        printf("supports only png screenshots");
        return;
    }

    auto wh = yglu::int2{0, 0};
    auto pixels = yglu::ui::get_screenshot(win, wh);
    yimg::save_image(imfilename, wh[0], wh[1], 4, nullptr,
                     (unsigned char*)pixels.data());
}

void draw_scene(yglu::ui::window* win) {
    auto st = (state*)yglu::ui::get_user_pointer(win);
    auto pars = st->pars;
    yapp::shade_scene(st->scene, pars, st->shst);
}

// deprecated function removed from build
void draw_debug(yglu::ui::window* win) {
#if 0
    auto pars = (yisym_app::params*)glfwGetWindowUserPointer(window);

    if (!pars->show_debug) return;

    // begin frame
    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    auto cam = pars->scene->cameras[pars->camera_id];
    auto camera_xform = ym::to_mat((ym::frame3f)cam->frame);
    auto camera_view = ym::to_mat(ym::inverse((ym::frame3f)cam->frame));
    auto camera_proj =
        ym::perspective_mat4(cam->yfov, cam->aspect, 0.1f, 10000.0f);

    yglu::legacy::begin_frame(*(yglu::float4x4*)&camera_xform,
                              *(yglu::float4x4*)&camera_view,
                              *(yglu::float4x4*)&camera_proj, false, true);

    for (auto i = 0; i < pars->scene->shapes.size(); i++) {
        auto shape = pars->scene->shapes[i];
        // auto& body = rigid_scene.shapes[i];

        auto xform = ym::to_mat((ym::frame3f)shape->frame);
        yglu::legacy::begin_shape(*(yglu::float4x4*)&xform);

        yglu::legacy::set_material({0, 0, 0}, {0, 0, 0}, {0, 0, 0}, 1, 0, true);

        glLineWidth(2);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        if (!shape->points.empty()) {
            yglu::legacy::draw_elems((int)shape->points.size(),
                                     shape->points.data(), yglu::etype::point,
                                     shape->pos.data(), nullptr, nullptr,
                                     nullptr);
        } else if (!shape->lines.empty()) {
            yglu::legacy::draw_elems((int)shape->lines.size(),
                                     (int*)shape->lines.data(),
                                     yglu::etype::line, shape->pos.data(),
                                     nullptr, nullptr, nullptr);
        } else if (!shape->triangles.empty()) {
            yglu::legacy::draw_elems((int)shape->triangles.size(),
                                     (int*)shape->triangles.data(),
                                     yglu::etype::triangle, shape->pos.data(),
                                     nullptr, nullptr, nullptr);
        }

        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glLineWidth(1);

        yglu::legacy::end_shape();
    }

    for (auto collision : pars->rigid_scene->__collisions) {
        int point[] = {0}, line[] = {0, 1};
        yglu::legacy::begin_shape(*(yglu::float4x4*)&ym::identity_mat4f);

        yglu::legacy::set_material({1, 0, 0}, {0, 0, 0}, {0, 0, 0}, 1, 0,
                                   false);

        ym::vec3f pos[] = {ym::pos(collision.frame),
                           ym::pos(collision.frame) +
                               collision.depth * collision.frame[2]};

        glPointSize(10);
        glLineWidth(4);
        yglu::legacy::draw_elems(1, point, yglu::etype::point,
                                 (yglu::float3*)pos, nullptr, nullptr, nullptr);
        yglu::legacy::draw_elems(1, line, yglu::etype::line, (yglu::float3*)pos,
                                 nullptr, nullptr, nullptr);
        glLineWidth(1);
        glPointSize(1);

        yglu::legacy::set_material({0, 1, 0}, {0, 0, 0}, {0, 0, 0}, 1, 0,
                                   false);

        ym::vec3f posi[] = {ym::pos(collision.frame),
                            ym::pos(collision.frame) +
                                collision.impulse * 10.0f};
        glLineWidth(4);
        yglu::legacy::draw_elems(1, line, yglu::etype::line,
                                 (yglu::float3*)posi, nullptr, nullptr,
                                 nullptr);
        glLineWidth(1);

        yglu::legacy::set_material({0.5, 0.5, 1}, {0, 0, 0}, {0, 0, 0}, 1, 0,
                                   false);

        ym::vec3f posj[] = {ym::pos(collision.frame),
                            ym::pos(collision.frame) +
                                pars->dt * collision.vel_before};
        glLineWidth(4);
        yglu::legacy::draw_elems(1, line, yglu::etype::line,
                                 (yglu::float3*)posj, nullptr, nullptr,
                                 nullptr);
        glLineWidth(1);

        yglu::legacy::set_material({0, 1, 1}, {0, 0, 0}, {0, 0, 0}, 1, 0,
                                   false);

        ym::vec3f posk[] = {ym::pos(collision.frame),
                            ym::pos(collision.frame) +
                                pars->dt * collision.vel_after};
        glLineWidth(4);
        yglu::legacy::draw_elems(1, line, yglu::etype::line,
                                 (yglu::float3*)posk, nullptr, nullptr,
                                 nullptr);
        glLineWidth(1);

        yglu::legacy::end_shape();
    }

    yglu::legacy::end_frame();

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glEnable(GL_DEPTH_TEST);
#endif

#if 0
    static int point[] = {0}, line[] = {0, 1};
    for (int i = 0; i < rigid_scene->ncollisions; i++) {
        ysym::_collision* col = rigid_scene->collisions + i;
        yocto::ym::mat4f identity4f = yocto::ym::identity_mat4f;
        yglu::stdshader_begin_shape(prog, identity4f.data());

        float zero[3] = {0, 0, 0};
        float red[3] = {1, 0, 0}, green[3] = {0, 1, 0}, cyan[3] = {0, 1, 1},
              yellow[3] = {1, 1, 0};

        yglu::stdshader_set_material(prog, 1, red, zero, zero, 0, 0, 0, 0, 0,
                                     false);

        yocto::ym::vec3f pos[] = {
            col->frame.pos, col->frame.pos + col->depth * col->frame.norm};
        glPointSize(10);
        glLineWidth(4);
        yglu::stdshader_set_vert(prog, &pos->x, 0, 0, 0);
        yglu::stdshader_draw_elem(prog, 1, point, 1);
        yglu::stdshader_draw_elem(prog, 1, line, 2);
        glLineWidth(1);
        glPointSize(1);

        yglu::stdshader_set_material(prog, 1, green, zero, zero, 0, 0, 0, 0, 0,
                                     false);

        yocto::ym::vec3f posi[] = {col->frame.pos,
                                   col->frame.pos + 25.0f * col->impulse};
        glLineWidth(4);
        yglu::stdshader_set_vert(prog, &posi->x, 0, 0, 0);
        yglu::stdshader_draw_elem(prog, 1, line, 2);
        glLineWidth(1);

        yglu::stdshader_set_material(prog, 1, cyan, zero, zero, 0, 0, 0, 0, 0,
                                     false);

        yocto::ym::vec3f posj[] = {col->frame.pos,
                                   col->frame.pos + dt * col->vel_before};
        glLineWidth(4);
        yglu::stdshader_set_vert(prog, &posj->x, 0, 0, 0);
        yglu::stdshader_draw_elem(prog, 1, line, 2);
        glLineWidth(1);

        yglu::stdshader_set_material(prog, 1, yellow, zero, zero, 0, 0, 0, 0, 0,
                                     false);

        yocto::ym::vec3f posk[] = {col->frame.pos,
                                   col->frame.pos + dt * col->vel_after};
        glLineWidth(4);
        yglu::stdshader_set_vert(prog, &posk->x, 0, 0, 0);
        yglu::stdshader_draw_elem(prog, 1, line, 2);
        glLineWidth(1);

        yglu::stdshader_end_shape();
    }

#endif
}

void text_callback(yglu::ui::window* win, unsigned int key) {
    auto st = (state*)yglu::ui::get_user_pointer(win);
    auto pars = st->pars;
    switch (key) {
        case ' ': st->simulating = !st->simulating; break;
        case '/': {
            for (int sid = 0; sid < st->scene->shapes.size(); sid++) {
                st->scene->shapes[sid]->frame = st->initial_state[sid];
                ysym::set_body_frame(st->rigid_scene, sid,
                                     st->initial_state[sid]);
            }
            st->frame = 0;
        } break;
        case '.':
            yapp::simulate_step(st->scene, st->rigid_scene, pars->dt);
            st->frame += 1;
            break;
        case '[': pars->exposure -= 1; break;
        case ']': pars->exposure += 1; break;
        case '{': pars->gamma -= 0.1f; break;
        case '}': pars->gamma += 0.1f; break;
        case 'e': pars->edges = !pars->edges; break;
        case 'w': pars->wireframe = !pars->wireframe; break;
        case 's': save_screenshot(win, pars->imfilename); break;
        case 'C':
            pars->render_params.camera_id =
                (pars->render_params.camera_id + 1) % st->scene->cameras.size();
            break;
        case 'd': st->show_debug = !st->show_debug; break;
        default: printf("unsupported key\n"); break;
    }
}

void draw_widgets(yglu::ui::window* win) {
    auto st = (state*)yglu::ui::get_user_pointer(win);
    auto pars = st->pars;
    if (yglu::ui::begin_widgets(win)) {
        yglu::ui::dynamic_widget_layout(win, 1);
        yglu::ui::label_widget(win, pars->filenames[0]);
        yglu::ui::dynamic_widget_layout(win, 2);
        yglu::ui::int_label_widget(win, "frame", st->frame);
        yglu::ui::float_widget(win, "dt", &pars->dt, 0, 1, 1 / 240.0f);
        if (yglu::ui::button_widget(win, "start")) st->simulating = true;
        if (yglu::ui::button_widget(win, "stop")) st->simulating = false;
        if (yglu::ui::button_widget(win, "step")) {
            yapp::simulate_step(st->scene, st->rigid_scene, pars->dt);
            st->frame += 1;
        }
        if (yglu::ui::button_widget(win, "reset")) {
            for (int sid = 0; sid < st->scene->shapes.size(); sid++) {
                st->scene->shapes[sid]->frame = st->initial_state[sid];
                ysym::set_body_frame(st->rigid_scene, sid,
                                     st->initial_state[sid]);
                ybvh::set_shape_frame(st->scene_bvh, sid,
                                      st->initial_state[sid]);
            }
            st->frame = 0;
        }
        yglu::ui::int_widget(win, "camera", &pars->render_params.camera_id, 0,
                             (int)st->scene->cameras.size() - 1, 1);
        yglu::ui::bool_widget(win, "wireframe", &pars->wireframe);
        yglu::ui::bool_widget(win, "edges", &pars->edges);
        yglu::ui::dynamic_widget_layout(win, 1);
        yglu::ui::float_widget(win, "exposure", &pars->exposure, -20, 20, 1);
        yglu::ui::float_widget(win, "gamma", &pars->gamma, 0.1, 5, 0.1);
        if (yglu::ui::button_widget(win, "tesselate")) {
            for (auto shape : st->scene->shapes) {
                yshape::tesselate_stdshape(
                    (std::vector<yshape::int2>&)shape->lines,
                    (std::vector<yshape::int3>&)shape->triangles,
                    (std::vector<yshape::float3>&)shape->pos,
                    (std::vector<yshape::float3>&)shape->norm,
                    (std::vector<yshape::float2>&)shape->texcoord,
                    (std::vector<yshape::float3>&)shape->color, shape->radius);
            }
        }
    }
    yglu::ui::end_widgets(win);
}

void window_refresh_callback(yglu::ui::window* win) {
    draw_scene(win);
    draw_debug(win);
    draw_widgets(win);
    yglu::ui::swap_buffers(win);
}

void run_ui(state* st) {
    auto pars = st->pars;

    // window
    auto win = yglu::ui::init_window(pars->width, pars->height, "ysym",
                                     pars->legacy_gl, st);
    yglu::ui::set_callbacks(win, text_callback, window_refresh_callback);

    // window values
    int mouse_button = 0;
    ym::vec2f mouse_pos, mouse_last;

    // load textures
    st->shst = yapp::init_shade_state(st->scene, pars);

    yglu::ui::init_widgets(win);

    while (!yglu::ui::should_close(win)) {
        mouse_last = mouse_pos;
        mouse_pos = yglu::ui::get_mouse_posf(win);
        mouse_button = yglu::ui::get_mouse_button(win);

        yglu::ui::set_window_title(win, ("yshade | " + pars->filenames[0]));

        // handle mouse
        if (mouse_button && mouse_pos != mouse_last &&
            !yglu::ui::get_widget_active(win)) {
            auto dolly = 0.0f;
            auto pan = ym::zero2f;
            auto rotate = ym::zero2f;
            switch (mouse_button) {
                case 1: rotate = (mouse_pos - mouse_last) / 100.0f; break;
                case 2: dolly = (mouse_pos[0] - mouse_last[0]) / 100.0f; break;
                case 3: pan = (mouse_pos - mouse_last) / 100.0f; break;
                default: break;
            }

            auto cam = st->scene->cameras[pars->render_params.camera_id];
            ym::turntable((ym::frame3f&)cam->frame, cam->focus, rotate, dolly,
                          pan);
        }

        // draw
        draw_scene(win);
        draw_debug(win);
        draw_widgets(win);

        // swap buffers
        yglu::ui::swap_buffers(win);

        // advance if simulating
        if (st->simulating) {
            yapp::simulate_step(st->scene, st->rigid_scene, pars->dt);
            st->frame += 1;
        }

        // check for screenshot
        if (pars->no_ui) {
            save_screenshot(win, pars->imfilename);
            break;
        }

        // event hadling
        if (st->simulating)
            yglu::ui::poll_events(win);
        else
            yglu::ui::wait_events(win);
    }

    yglu::ui::clear_widgets(win);
    yglu::ui::clear_window(win);
}

int main(int argc, char* argv[]) {
    // params
    auto pars = yapp::init_params("interactively simulate scenes", argc, argv,
                                  false, true, true, true);

    // init state
    auto st = new state();
    st->pars = pars;

    // setting up rendering
    st->scene = yapp::load_scenes(pars->filenames, pars->scene_scale);
    st->scene_bvh = yapp::make_bvh(st->scene);
    st->rigid_scene = yapp::make_rigid_scene(st->scene, st->scene_bvh);

    // initialize simulation
    ysym::init_simulation(st->rigid_scene);

    // save init values
    st->initial_state.resize(st->scene->shapes.size());
    for (auto i = 0; i < st->initial_state.size(); i++)
        st->initial_state[i] = st->scene->shapes[i]->frame;

    // run ui
    run_ui(st);

    // cleanup
    delete st;

    // done
    return 0;
}
