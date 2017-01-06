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

#include "ytrace.h"
#include "yui.h"

#include "../yocto/yocto_glu.h"

namespace yitrace_app {

// interative params
struct params : virtual ytrace_app::params {
    // rendered image
    std::vector<ym::vec4b> ldr;

    // gl
    bool no_ui = false;
    bool legacy_gl = false;

    // progressive rendering
    int preview_width = 0;
    int preview_height = 0;
    std::vector<ym::vec4f> preview;
    int cur_sample = 0, cur_block = 0;
    int blocks_per_update = 8;

    // view variables
    float background = 0;

    // shading
    bool scene_updated = true;
    yglu::uint texture_id = 0;
    float texture_exposure, texture_gamma;

    // widgets
    void* widget_ctx = nullptr;
};

// init params
void init_params(yitrace_app::params* pars, ycmd::parser* parser) {
    // parse cmdline
    pars->no_ui =
        ycmd::parse_flag(parser, "--no-ui", "", "runs offline", false);
    pars->legacy_gl = ycmd::parse_flag(parser, "--legacy_opengl", "-L",
                                       "uses legacy OpenGL", false);

    // init params
    ytrace_app::init_params(pars, parser);

    // image rendering params
    pars->ldr.resize(pars->width * pars->height, {0, 0, 0, 0});
    pars->preview_width = pars->width / pars->block_size;
    pars->preview_height = pars->height / pars->block_size;
    pars->preview.resize(pars->preview_width * pars->preview_height,
                         {0, 0, 0, 0});
    pars->texture_exposure = pars->exposure;
    pars->texture_gamma = pars->gamma;
}
};

// nuklear
const int hud_width = 256;

void text_callback(GLFWwindow* window, unsigned int key) {
    auto pars = (yitrace_app::params*)glfwGetWindowUserPointer(window);
    auto nuklear_ctx = (nk_context*)pars->widget_ctx;

    nk_glfw3_gl3_char_callback(window, key);
    if (nk_item_is_any_active(nuklear_ctx)) return;
    switch (key) {
        case '[': pars->exposure -= 1; break;
        case ']': pars->exposure += 1; break;
        case '{': pars->gamma -= 0.1f; break;
        case '}': pars->gamma += 0.1f; break;
        case '1':
            pars->exposure = 0;
            pars->gamma = 1;
            break;
        case '2':
            pars->exposure = 0;
            pars->gamma = 2.2f;
            break;
        case 's':
            ytrace_app::save_image(pars->imfilename, pars->width, pars->height,
                                   pars->hdr.data(), pars->exposure,
                                   pars->gamma);
            break;
        default: printf("unsupported key\n"); break;
    }
}

void draw_image(GLFWwindow* window) {
    auto pars = (yitrace_app::params*)glfwGetWindowUserPointer(window);
    auto framebuffer_size = yui::framebuffer_size(window);
    glViewport(0, 0, framebuffer_size[0], framebuffer_size[1]);

    // begin frame
    glClearColor(pars->background, pars->background, pars->background, 0);
    glDisable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT);

    // draw image
    auto window_size = yui::window_size(window);
    if (pars->legacy_gl) {
        yglu::legacy::draw_image(pars->texture_id, pars->width, pars->height,
                                 window_size[0], window_size[1], 0, 0, 1);
    } else {
        yglu::modern::shade_image(pars->texture_id, pars->width, pars->height,
                                  window_size[0], window_size[1], 0, 0, 1);
    }
}

void draw_widgets(GLFWwindow* window) {
    auto pars = (yitrace_app::params*)glfwGetWindowUserPointer(window);
    auto nuklear_ctx = (nk_context*)pars->widget_ctx;
    auto window_size = yui::window_size(window);
    if (pars->legacy_gl) {
        nk_glfw3_gl2_new_frame();
    } else {
        nk_glfw3_gl3_new_frame();
    }
    if (nk_begin(nuklear_ctx, "ytrace", nk_rect(window_size[0] - hud_width, 0,
                                                hud_width, window_size[1]),
                 NK_WINDOW_BORDER | NK_WINDOW_MOVABLE | NK_WINDOW_SCALABLE |
                     NK_WINDOW_MINIMIZABLE | NK_WINDOW_TITLE)) {
        nk_layout_row_dynamic(nuklear_ctx, 30, 1);
        nk_label(nuklear_ctx, pars->filename.c_str(), NK_TEXT_LEFT);
        nk_layout_row_dynamic(nuklear_ctx, 30, 3);
        nk_value_int(nuklear_ctx, "w", pars->width);
        nk_value_int(nuklear_ctx, "h", pars->height);
        nk_value_int(nuklear_ctx, "s", pars->cur_sample);
        nk_layout_row_dynamic(nuklear_ctx, 30, 1);
        nk_property_int(nuklear_ctx, "samples", 0,
                        &pars->render_params.nsamples, 1000000, 1, 1);
        nk_layout_row_dynamic(nuklear_ctx, 30, 2);
        nk_property_int(nuklear_ctx, "camera", 0,
                        &pars->render_params.camera_id,
                        (int)pars->scene->cameras.size() - 1, 1, 1);
        nk_layout_row_dynamic(nuklear_ctx, 30, 1);
        nk_property_float(nuklear_ctx, "exposure", -20, &pars->exposure, 20, 1,
                          1);
        nk_property_float(nuklear_ctx, "gamma", 0.1, &pars->gamma, 5, 0.1, 0.1);
    }
    nk_end(nuklear_ctx);

    if (pars->legacy_gl) {
        nk_glfw3_gl2_render(NK_ANTI_ALIASING_ON, 512 * 1024, 128 * 1024);
    } else {
        nk_glfw3_gl3_render(NK_ANTI_ALIASING_ON, 512 * 1024, 128 * 1024);
    }
}

void window_refresh_callback(GLFWwindow* window) {
    draw_image(window);
    draw_widgets(window);
    glfwSwapBuffers(window);
}

bool update(yitrace_app::params* pars) {
    if (pars->scene_updated) {
        // update camera
        auto cam = pars->scene->cameras[pars->render_params.camera_id];
        auto trace_cam =
            pars->trace_scene->cameras[pars->render_params.camera_id];
        *trace_cam = {cam->frame, cam->yfov, cam->aspect, cam->aperture,
                      cam->focus};

        // render preview
        auto pparams = pars->render_params;
        pparams.nsamples = 1;
        ytrace::trace_image(pars->trace_scene, pars->preview_width,
                            pars->preview_height,
                            (ytrace::float4*)pars->preview.data(), pparams);
        for (auto qj = 0; qj < pars->preview_height; qj++) {
            for (auto qi = 0; qi < pars->preview_width; qi++) {
                for (auto j = qj * pars->block_size;
                     j < ym::min((qj + 1) * pars->block_size, pars->height);
                     j++) {
                    for (auto i = qi * pars->block_size;
                         i < ym::min((qi + 1) * pars->block_size, pars->width);
                         i++) {
                        pars->hdr[j * pars->width + i] =
                            pars->preview[qj * pars->preview_width + qi];
                    }
                }
            }
        }
        ym::exposure_gamma(
            pars->width, pars->height, 4, (const float*)pars->hdr.data(),
            (unsigned char*)pars->ldr.data(), pars->exposure, pars->gamma);
        if (pars->legacy_gl) {
            yglu::legacy::update_texture(
                pars->texture_id, pars->width, pars->height, 4,
                (unsigned char*)pars->ldr.data(), false);
        } else {
            yglu::modern::update_texture(
                pars->texture_id, pars->width, pars->height, 4,
                (unsigned char*)pars->ldr.data(), false);
        }

        // reset current counters
        pars->cur_sample = 0;
        pars->cur_block = 0;
        pars->scene_updated = false;
    } else {
        if (pars->cur_sample == pars->render_params.nsamples) return false;
        std::vector<std::future<void>> futures;
        for (auto b = 0; pars->cur_block < pars->blocks.size() &&
                         b < pars->blocks_per_update;
             pars->cur_block++, b++) {
            auto block = pars->blocks[pars->cur_block];
            futures.push_back(pars->pool->enqueue([block, pars]() {
                ytrace::trace_block(
                    pars->trace_scene, pars->width, pars->height,
                    (ytrace::float4*)pars->hdr.data(), block[0], block[1],
                    block[2], block[3], pars->cur_sample, pars->cur_sample + 1,
                    pars->render_params, true);
                ym::exposure_gamma(pars->width, pars->height, 4,
                                   (const float*)pars->hdr.data(),
                                   (unsigned char*)pars->ldr.data(),
                                   pars->exposure, pars->gamma, block[0],
                                   block[1], block[2], block[3]);
            }));
        }
        for (auto& future : futures) future.wait();
        if (pars->texture_exposure != pars->exposure ||
            pars->texture_gamma != pars->gamma) {
            ym::exposure_gamma(
                pars->width, pars->height, 4, (const float*)pars->hdr.data(),
                (unsigned char*)pars->ldr.data(), pars->exposure, pars->gamma);
            pars->texture_exposure = pars->exposure;
            pars->texture_gamma = pars->gamma;
        }
        if (pars->legacy_gl) {
            yglu::legacy::update_texture(
                pars->texture_id, pars->width, pars->height, 4,
                (unsigned char*)pars->ldr.data(), false);
        } else {
            yglu::modern::update_texture(
                pars->texture_id, pars->width, pars->height, 4,
                (unsigned char*)pars->ldr.data(), false);
        }
        if (pars->cur_block == pars->blocks.size()) {
            pars->cur_block = 0;
            pars->cur_sample++;
        }
    }

    return true;
}

void run_ui(yitrace_app::params* pars) {
    // window
    auto window = yui::init_glfw(pars->width, pars->height, "ytrace",
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

    pars->widget_ctx = yui::init_nuklear(window, pars->legacy_gl);

    if (pars->legacy_gl) {
        pars->texture_id = yglu::legacy::make_texture(
            pars->width, pars->height, 4, (unsigned char*)pars->ldr.data(),
            false, false);
    } else {
        pars->texture_id = yglu::modern::make_texture(
            pars->width, pars->height, 4, (unsigned char*)pars->ldr.data(),
            false, false, false);
    }

    while (!glfwWindowShouldClose(window)) {
        glfwGetWindowSize(window, &window_size[0], &window_size[1]);
        glfwGetFramebufferSize(window, &framebuffer_size[0],
                               &framebuffer_size[1]);

        mouse_last = mouse_pos;
        mouse_pos = yui::mouse_pos(window);
        mouse_button = yui::mouse_button(window);

        glfwSetWindowTitle(window, ("ytrace | " + pars->filename + " | " +
                                    std::to_string(pars->width) + "x" +
                                    std::to_string(pars->height) + "@" +
                                    std::to_string(pars->cur_sample))
                                       .c_str());

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

            auto cam = pars->scene->cameras[pars->render_params.camera_id];
            ym::turntable((ym::frame3f&)cam->frame, cam->focus, rotate, dolly,
                          pan);
            pars->scene_updated = true;
        }

        // update
        auto updated = update(pars);

        // draw
        draw_image(window);

        // make ui
        draw_widgets(window);

        // swap buffers
        glfwSwapBuffers(window);

        // event hadling
        if (updated)
            glfwPollEvents();
        else
            glfwWaitEvents();
    }

    yui::clear_nuklear((nk_context*)pars->widget_ctx, pars->legacy_gl);
    yui::clear_glfw(window);
}

int main(int argc, char* argv[]) {
    // params
    auto pars = new yitrace_app::params();
    auto parser = ycmd::make_parser(argc, argv, "trace meshes");
    yitrace_app::init_params(pars, parser);
    ycmd::check_parser(parser);

    // launching renderer
    if (pars->no_ui) {
        ytrace_app::render(pars);
    } else {
        // run ui
        run_ui(pars);
    }

    // done
    delete pars;
    return EXIT_SUCCESS;
}
