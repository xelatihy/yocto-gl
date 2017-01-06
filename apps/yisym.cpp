//
// LICENSE:
//
// Copyright (c) 2016 Fabio Pellacini
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

#include "yshade.h"
#include "ysym.h"
#include "yui.h"

namespace yisym_app {

struct params : virtual ysym_app::params, virtual yshade_app::params {
    bool simulating = false;
    int frame = 0;
    std::vector<ym::frame3f> initial_state;
    bool show_debug = false;
};

void init_params(params* pars, ycmd::parser* parser) {
    // init params
    ysym_app::init_params(pars, parser);
    yshade_app::init_params(pars, parser);

    // init values
    // save out init state
    pars->initial_state.resize(pars->scene->shapes.size());
    for (auto i = 0; i < pars->initial_state.size(); i++)
        pars->initial_state[i] = pars->scene->shapes[i]->frame;
}
}  // namespace

const int hud_width = 256;

void save_screenshot(GLFWwindow* window, const std::string& imfilename) {
    if (ycmd::get_extension(imfilename) != ".png") {
        printf("supports only png screenshots");
        return;
    }

    auto wh = yui::framebuffer_size(window);
    auto pixels = std::vector<unsigned char>(wh[0] * wh[1] * 4);
    glReadBuffer(GL_FRONT);
    glReadPixels(0, 0, wh[0], wh[1], GL_RGBA, GL_UNSIGNED_BYTE, pixels.data());
    std::vector<unsigned char> line(wh[0] * 4);
    for (int j = 0; j < wh[1] / 2; j++) {
        memcpy(line.data(), pixels.data() + j * wh[0] * 4, wh[0] * 4);
        memcpy(pixels.data() + j * wh[0] * 4,
               pixels.data() + (wh[1] - 1 - j) * wh[0] * 4, wh[0] * 4);
        memcpy(pixels.data() + (wh[1] - 1 - j) * wh[0] * 4, line.data(),
               wh[0] * 4);
    }
    stbi_write_png(imfilename.c_str(), wh[0], wh[1], 4, pixels.data(),
                   wh[0] * 4);
}

void draw_scene(GLFWwindow* window) {
    auto pars = (yisym_app::params*)glfwGetWindowUserPointer(window);
    yshade_app::render(pars);
}

// deprecated function removed from build
void draw_debug(GLFWwindow* window) {
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

void text_callback(GLFWwindow* window, unsigned int key) {
    auto pars = (yisym_app::params*)glfwGetWindowUserPointer(window);
    auto nuklear_ctx = (nk_context*)pars->widget_ctx;

    nk_glfw3_gl3_char_callback(window, key);
    if (nk_item_is_any_active(nuklear_ctx)) return;
    switch (key) {
        case ' ': pars->simulating = !pars->simulating; break;
        case '/': {
            for (int sid = 0; sid < pars->scene->shapes.size(); sid++) {
                pars->scene->shapes[sid]->frame = pars->initial_state[sid];
                ysym::set_body_frame(pars->rigid_scene, sid,
                                     pars->initial_state[sid]);
            }
            pars->frame = 0;
        } break;
        case '.':
            ysym_app::simulate_step(pars->scene, pars->rigid_scene, pars->dt);
            pars->frame += 1;
            break;
        case '[': pars->exposure -= 1; break;
        case ']': pars->exposure += 1; break;
        case '{': pars->gamma -= 0.1f; break;
        case '}': pars->gamma += 0.1f; break;
        case 'e': pars->edges = !pars->edges; break;
        case 'w': pars->wireframe = !pars->wireframe; break;
        case 's': save_screenshot(window, pars->imfilename); break;
        case 'c': pars->camera_lights = !pars->camera_lights; break;
        case 'C':
            pars->camera_id =
                (pars->camera_id + 1) % pars->scene->cameras.size();
            break;
        case 'd': pars->show_debug = !pars->show_debug; break;
        default: printf("unsupported key\n"); break;
    }
}

void draw_widgets(GLFWwindow* window) {
    auto pars = (yisym_app::params*)glfwGetWindowUserPointer(window);
    auto nuklear_ctx = (nk_context*)pars->widget_ctx;
    auto window_size = yui::window_size(window);
    if (pars->legacy_gl) {
        nk_glfw3_gl2_new_frame();
    } else {
        nk_glfw3_gl3_new_frame();
    }
    if (nk_begin(nuklear_ctx, "ysym", nk_rect(window_size[0] - hud_width, 0,
                                              hud_width, window_size[1]),
                 NK_WINDOW_BORDER | NK_WINDOW_MOVABLE | NK_WINDOW_SCALABLE |
                     NK_WINDOW_MINIMIZABLE | NK_WINDOW_TITLE)) {
        nk_layout_row_dynamic(nuklear_ctx, 30, 1);
        nk_label(nuklear_ctx, pars->filename.c_str(), NK_TEXT_LEFT);
        nk_layout_row_dynamic(nuklear_ctx, 30, 2);
        nk_value_int(nuklear_ctx, "frame", pars->frame);
        nk_property_float(nuklear_ctx, "dt", 0, &pars->dt, 1, 1 / 240.0f,
                          1 / 240.0f);
        if (nk_button_label(nuklear_ctx, "start")) pars->simulating = true;
        if (nk_button_label(nuklear_ctx, "stop")) pars->simulating = false;
        if (nk_button_label(nuklear_ctx, "step")) {
            ysym_app::simulate_step(pars->scene, pars->rigid_scene, pars->dt);
            pars->frame += 1;
        }
        if (nk_button_label(nuklear_ctx, "reset")) {
            for (int sid = 0; sid < pars->scene->shapes.size(); sid++) {
                pars->scene->shapes[sid]->frame = pars->initial_state[sid];
                pars->rigid_scene->shapes[sid]->frame =
                    pars->initial_state[sid];
            }
            pars->frame = 0;
        }
        nk_property_int(nuklear_ctx, "camera", 0, &pars->camera_id,
                        (int)pars->scene->cameras.size() - 1, 1, 1);
        pars->camera_lights =
            nk_check_label(nuklear_ctx, "eyelight", pars->camera_lights);
        pars->wireframe =
            nk_check_label(nuklear_ctx, "wireframe", pars->wireframe);
        pars->edges = nk_check_label(nuklear_ctx, "edges", pars->edges);
        nk_layout_row_dynamic(nuklear_ctx, 30, 1);
        nk_property_float(nuklear_ctx, "exposure", -20, &pars->exposure, 20, 1,
                          1);
        nk_property_float(nuklear_ctx, "gamma", 0.1, &pars->gamma, 5, 0.1, 0.1);
        if (nk_button_label(nuklear_ctx, "tesselate")) {
            for (auto shape : pars->scene->shapes) {
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
    nk_end(nuklear_ctx);

    if (pars->legacy_gl) {
        nk_glfw3_gl2_render(NK_ANTI_ALIASING_ON, 512 * 1024, 128 * 1024);
    } else {
        nk_glfw3_gl3_render(NK_ANTI_ALIASING_ON, 512 * 1024, 128 * 1024);
    }
}

void window_refresh_callback(GLFWwindow* window) {
    draw_scene(window);
    draw_debug(window);
    draw_widgets(window);
    glfwSwapBuffers(window);
}

void run_ui(yisym_app::params* pars) {
    // window
    auto window = yui::init_glfw(pars->width, pars->height, "ysym",
                                 pars->legacy_gl, pars, text_callback);

    // callbacks
    glfwSetWindowRefreshCallback(window, window_refresh_callback);
    glfwSetScrollCallback(window, nk_gflw3_scroll_callback);

    // window values
    int mouse_button = 0;
    ym::vec2f mouse_pos, mouse_last;
    ym::vec2i window_size, framebuffer_size;

// init gl extensions
#ifndef __APPLE__
    if (glewInit() != GLEW_OK) exit(EXIT_FAILURE);
#endif

    // load textures
    yshade_app::init(pars);

    pars->widget_ctx = yui::init_nuklear(window, pars->legacy_gl);

    while (!glfwWindowShouldClose(window)) {
        glfwGetWindowSize(window, &window_size[0], &window_size[1]);
        glfwGetFramebufferSize(window, &framebuffer_size[0],
                               &framebuffer_size[1]);

        mouse_last = mouse_pos;
        mouse_pos = yui::mouse_pos(window);
        mouse_button = yui::mouse_button(window);

        glfwSetWindowTitle(window, ("yshade | " + pars->filename).c_str());

        // handle mouse
        if (mouse_button && mouse_pos != mouse_last &&
            !nk_item_is_any_active((nk_context*)pars->widget_ctx)) {
            auto dolly = 0.0f;
            auto pan = ym::zero2f;
            auto rotate = ym::zero2f;
            switch (mouse_button) {
                case 1: rotate = (mouse_pos - mouse_last) / 100.0f; break;
                case 2: dolly = (mouse_pos[0] - mouse_last[0]) / 100.0f; break;
                case 3: pan = (mouse_pos - mouse_last) / 100.0f; break;
                default: break;
            }

            auto cam = pars->scene->cameras[pars->camera_id];
            ym::turntable((ym::frame3f&)cam->frame, cam->focus, rotate, dolly,
                          pan);
        }

        // draw
        draw_scene(window);
        draw_debug(window);
        draw_widgets(window);

        // swap buffers
        glfwSwapBuffers(window);

        // advance if simulating
        if (pars->simulating) {
            ysym_app::simulate_step(pars->scene, pars->rigid_scene, pars->dt);
            pars->frame += 1;
        }

        // check for screenshot
        if (pars->no_ui) {
            save_screenshot(window, pars->imfilename);
            break;
        }

        // event hadling
        if (pars->simulating)
            glfwPollEvents();
        else
            glfwWaitEvents();
    }

    yui::clear_nuklear((nk_context*)pars->widget_ctx, pars->legacy_gl);
    yui::clear_glfw(window);
}

int main(int argc, char* argv[]) {
    // command line
    auto pars = new yisym_app::params();
    auto parser = ycmd::make_parser(argc, argv, "rigid body simulation");
    yisym_app::init_params(pars, parser);
    ycmd::check_parser(parser);

    // run ui
    run_ui(pars);

    // done
    delete pars;
    return 0;
}
