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
#include "yui.h"

#include "../yocto/yocto_glu.h"

#include "ThreadPool.h"

// interative state
struct state {
    // params
    yapp::params* pars = nullptr;

    // scene
    yapp::scene* scene = nullptr;
    ybvh::scene* scene_bvh = nullptr;
    ytrace::scene* trace_scene = nullptr;

    // rendered image
    std::vector<ym::vec4f> hdr;
    std::vector<ym::vec4b> ldr;

    // rendering aids
    ThreadPool* pool = nullptr;
    std::vector<yapp::int4> blocks;

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
    bool texture_srgb;

    // widgets
    void* widget_ctx = nullptr;
};

// nuklear
const int hud_width = 256;

void text_callback(GLFWwindow* window, unsigned int key) {
    auto st = (state*)glfwGetWindowUserPointer(window);
    auto pars = st->pars;
    auto nuklear_ctx = (nk_context*)st->widget_ctx;

    nk_glfw3_gl3_char_callback(window, key);
    if (nk_item_is_any_active(nuklear_ctx)) return;
    switch (key) {
        case '[': pars->exposure -= 1; break;
        case ']': pars->exposure += 1; break;
        case '{': pars->gamma -= 0.1f; break;
        case '}': pars->gamma += 0.1f; break;
        case '\\': pars->srgb = !pars->srgb; break;
        case '0':
            pars->exposure = 0;
            pars->gamma = 1;
            pars->srgb = false;
            break;
        case '1':
            pars->exposure = 0;
            pars->gamma = 1;
            pars->srgb = true;
            break;
        case '2':
            pars->exposure = 0;
            pars->gamma = 2.2f;
            pars->srgb = true;
            break;
        case 's':
            yapp::save_image(pars->imfilename, pars->width, pars->height,
                             st->hdr.data(), pars->exposure, pars->gamma,
                             pars->srgb);
            break;
        default: printf("unsupported key\n"); break;
    }
}

void draw_image(GLFWwindow* window) {
    auto st = (state*)glfwGetWindowUserPointer(window);
    auto pars = st->pars;
    auto framebuffer_size = yui::framebuffer_size(window);
    glViewport(0, 0, framebuffer_size[0], framebuffer_size[1]);

    // begin frame
    glClearColor(pars->background[0], pars->background[1], pars->background[2],
                 pars->background[3]);
    glDisable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT);

    // draw image
    auto window_size = yui::window_size(window);
    if (pars->legacy_gl) {
        yglu::legacy::draw_image(st->texture_id, pars->width, pars->height,
                                 window_size[0], window_size[1], 0, 0, 1);
    } else {
        yglu::modern::shade_image(st->texture_id, pars->width, pars->height,
                                  window_size[0], window_size[1], 0, 0, 1);
    }
}

void draw_widgets(GLFWwindow* window) {
    auto st = (state*)glfwGetWindowUserPointer(window);
    auto pars = st->pars;
    auto nuklear_ctx = (nk_context*)st->widget_ctx;
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
        nk_value_int(nuklear_ctx, "s", st->cur_sample);
        nk_layout_row_dynamic(nuklear_ctx, 30, 1);
        nk_property_int(nuklear_ctx, "samples", 0,
                        &pars->render_params.nsamples, 1000000, 1, 1);
        nk_layout_row_dynamic(nuklear_ctx, 30, 2);
        nk_property_int(nuklear_ctx, "camera", 0,
                        &pars->render_params.camera_id,
                        (int)st->scene->cameras.size() - 1, 1, 1);
        nk_layout_row_dynamic(nuklear_ctx, 30, 1);
        nk_property_float(nuklear_ctx, "hdr exposure", -20, &pars->exposure, 20,
                          1, 1);
        nk_property_float(nuklear_ctx, "hdr gamma", 0.1, &pars->gamma, 5, 0.1,
                          0.1);
        pars->srgb = nk_check_label(nuklear_ctx, "hdr srgb output", pars->srgb);
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

bool update(state* st) {
    auto pars = st->pars;
    if (st->scene_updated) {
        // update camera
        auto cam = st->scene->cameras[pars->render_params.camera_id];
        ytrace::set_camera(st->trace_scene, pars->render_params.camera_id,
                           cam->frame, cam->yfov, cam->aspect, cam->aperture,
                           cam->focus);

        // render preview
        auto pparams = pars->render_params;
        pparams.nsamples = 1;
        ytrace::trace_image(st->trace_scene, st->preview_width,
                            st->preview_height,
                            (ytrace::float4*)st->preview.data(), pparams);
        for (auto qj = 0; qj < st->preview_height; qj++) {
            for (auto qi = 0; qi < st->preview_width; qi++) {
                for (auto j = qj * pars->block_size;
                     j < ym::min((qj + 1) * pars->block_size, pars->height);
                     j++) {
                    for (auto i = qi * pars->block_size;
                         i < ym::min((qi + 1) * pars->block_size, pars->width);
                         i++) {
                        st->hdr[j * pars->width + i] =
                            st->preview[qj * st->preview_width + qi];
                    }
                }
            }
        }
        ym::exposure_gamma(pars->width, pars->height, 4,
                           (const float*)st->hdr.data(),
                           (unsigned char*)st->ldr.data(), pars->exposure,
                           pars->gamma, pars->srgb);
        if (pars->legacy_gl) {
            yglu::legacy::update_texture(st->texture_id, pars->width,
                                         pars->height, 4,
                                         (unsigned char*)st->ldr.data(), false);
        } else {
            yglu::modern::update_texture(st->texture_id, pars->width,
                                         pars->height, 4,
                                         (unsigned char*)st->ldr.data(), false);
        }

        // reset current counters
        st->cur_sample = 0;
        st->cur_block = 0;
        st->scene_updated = false;
    } else {
        if (st->cur_sample == pars->render_params.nsamples) return false;
        std::vector<std::future<void>> futures;
        for (auto b = 0;
             st->cur_block < st->blocks.size() && b < st->blocks_per_update;
             st->cur_block++, b++) {
            auto block = st->blocks[st->cur_block];
            futures.push_back(st->pool->enqueue([st, block, pars]() {
                ytrace::trace_block(st->trace_scene, pars->width, pars->height,
                                    (ytrace::float4*)st->hdr.data(), block[0],
                                    block[1], block[2], block[3],
                                    st->cur_sample, st->cur_sample + 1,
                                    pars->render_params);
                ym::exposure_gamma(
                    pars->width, pars->height, 4, (const float*)st->hdr.data(),
                    (unsigned char*)st->ldr.data(), pars->exposure, pars->gamma,
                    pars->srgb, block[0], block[1], block[2], block[3]);
            }));
        }
        for (auto& future : futures) future.wait();
        if (st->texture_exposure != pars->exposure ||
            st->texture_gamma != pars->gamma ||
            st->texture_srgb != pars->srgb) {
            ym::exposure_gamma(pars->width, pars->height, 4,
                               (const float*)st->hdr.data(),
                               (unsigned char*)st->ldr.data(), pars->exposure,
                               pars->gamma, pars->srgb);
            st->texture_exposure = pars->exposure;
            st->texture_gamma = pars->gamma;
            st->texture_srgb = pars->srgb;
        }
        if (pars->legacy_gl) {
            yglu::legacy::update_texture(st->texture_id, pars->width,
                                         pars->height, 4,
                                         (unsigned char*)st->ldr.data(), false);
        } else {
            yglu::modern::update_texture(st->texture_id, pars->width,
                                         pars->height, 4,
                                         (unsigned char*)st->ldr.data(), false);
        }
        if (st->cur_block == st->blocks.size()) {
            st->cur_block = 0;
            if (st->pars->save_progressive &&
                (st->cur_sample + 1) % pars->save_progressive == 0) {
                auto imfilename =
                    ycmd::get_dirname(pars->imfilename) +
                    ycmd::get_basename(pars->imfilename) +
                    ycmd::format_str(".%04d", st->cur_sample + 1) +
                    ycmd::get_extension(pars->imfilename);
                ycmd::log_msgf(ycmd::log_level_info, "ytrace",
                               "saving image %s", imfilename.c_str());
                yapp::save_image(imfilename, pars->width, pars->height,
                                 st->hdr.data(), pars->exposure, pars->gamma,
                                 pars->srgb);
            }
            st->cur_sample++;
        }
    }

    return true;
}

void run_ui(state* st) {
    auto pars = st->pars;

    // window
    auto window = yui::init_glfw(pars->width, pars->height, "ytrace",
                                 pars->legacy_gl, st, text_callback);

    // callbacks
    glfwSetWindowRefreshCallback(window, window_refresh_callback);
    glfwSetScrollCallback(window, nk_gflw3_scroll_callback);

    // window values
    int mouse_button = 0;
    ym::vec2f mouse_pos, mouse_last;
    ym::vec2i window_size, framebuffer_size;

    st->widget_ctx = yui::init_nuklear(window, pars->legacy_gl);

    if (pars->legacy_gl) {
        st->texture_id = yglu::legacy::make_texture(
            pars->width, pars->height, 4, (unsigned char*)st->ldr.data(), false,
            false);
    } else {
        st->texture_id = yglu::modern::make_texture(
            pars->width, pars->height, 4, (unsigned char*)st->ldr.data(), false,
            false, false);
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
                                    std::to_string(st->cur_sample))
                                       .c_str());

        // handle mouse
        if (mouse_button && mouse_pos != mouse_last &&
            !nk_item_is_any_active((nk_context*)st->widget_ctx)) {
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
            st->scene_updated = true;
        }

        // update
        auto updated = update(st);

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

    yui::clear_nuklear((nk_context*)st->widget_ctx, pars->legacy_gl);
    yui::clear_glfw(window);
}

int main(int argc, char* argv[]) {
    // params
    auto pars = yapp::init_params("render scene with path tracing", argc, argv,
                                  true, false, false, false);

    // init state
    auto st = new state();
    st->pars = pars;

    // setting up rendering
    st->scene = yapp::load_scene(pars->filename, pars->scene_scale);
    if (!pars->envmap_filename.empty())
        yapp::load_envmap(st->scene, pars->envmap_filename, pars->envmap_scale);
    st->scene_bvh = yapp::make_bvh(st->scene);
    st->trace_scene = yapp::make_trace_scene(st->scene, st->scene_bvh,
                                             pars->render_params.camera_id);

    // image rendering params
    st->hdr.resize(pars->width * pars->height, {0, 0, 0, 0});
    st->ldr.resize(pars->width * pars->height, {0, 0, 0, 0});
    st->preview_width = pars->width / pars->block_size;
    st->preview_height = pars->height / pars->block_size;
    st->preview.resize(st->preview_width * st->preview_height, {0, 0, 0, 0});
    st->texture_exposure = pars->exposure;
    st->texture_gamma = pars->gamma;
    st->texture_srgb = pars->srgb;
    st->scene_updated = true;

    // progressive rendering
    if (!pars->nthreads) pars->nthreads = std::thread::hardware_concurrency();
    st->blocks =
        yapp::make_trace_blocks(pars->width, pars->height, pars->block_size);
    st->pool = new ThreadPool(pars->nthreads);

    // fixing camera
    for (auto cam : st->scene->cameras)
        cam->aspect = (float)pars->width / (float)pars->height;

    // init renderer
    ytrace::init_lights(st->trace_scene);

    // run ui
    run_ui(st);

    // done
    delete st->scene;
    ybvh::free_scene(st->scene_bvh);
    ytrace::free_scene(st->trace_scene);
    delete st;
    delete pars;
    return 0;
}
