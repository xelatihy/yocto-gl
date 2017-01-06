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
#include "yui.h"

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
    auto pars = (yshade_app::params*)glfwGetWindowUserPointer(window);
    auto window_size = yui::window_size(window);
    pars->scene->cameras[pars->camera_id]->aspect =
        (float)window_size[0] / (float)window_size[1];
    yshade_app::render(pars);
}

void text_callback(GLFWwindow* window, unsigned int key) {
    auto pars = (yshade_app::params*)glfwGetWindowUserPointer(window);
    auto nuklear_ctx = (nk_context*)pars->widget_ctx;
    nk_glfw3_gl3_char_callback(window, key);
    if (nk_item_is_any_active(nuklear_ctx)) return;
    switch (key) {
        case '[': pars->exposure -= 1; break;
        case ']': pars->exposure += 1; break;
        case '{': pars->gamma -= 0.1f; break;
        case '}': pars->gamma += 0.1f; break;
        case '\\': pars->srgb = !pars->srgb; break;
        case 'w': pars->wireframe = !pars->wireframe; break;
        case 'e': pars->edges = !pars->edges; break;
        case 's': save_screenshot(window, pars->imfilename); break;
        case 'c': pars->camera_lights = !pars->camera_lights; break;
        case 'C':
            pars->camera_id =
                (pars->camera_id + 1) % pars->scene->cameras.size();
            break;
        case 't': {
            for (auto shape : pars->scene->shapes) {
                yshape::tesselate_stdshape(
                    (std::vector<yshape::int2>&)shape->lines,
                    (std::vector<yshape::int3>&)shape->triangles,
                    (std::vector<yshape::float3>&)shape->pos,
                    (std::vector<yshape::float3>&)shape->norm,
                    (std::vector<yshape::float2>&)shape->texcoord,
                    (std::vector<yshape::float3>&)shape->color, shape->radius);
            }
        } break;
        default: printf("unsupported key\n"); break;
    }
}

void draw_widgets(GLFWwindow* window) {
    auto pars = (yshade_app::params*)glfwGetWindowUserPointer(window);
    auto nuklear_ctx = (nk_context*)pars->widget_ctx;
    auto window_size = yui::window_size(window);
    if (pars->legacy_gl) {
        nk_glfw3_gl2_new_frame();
    } else {
        nk_glfw3_gl3_new_frame();
    }
    if (nk_begin(nuklear_ctx, "yshade", nk_rect(window_size[0] - hud_width, 0,
                                                hud_width, window_size[1]),
                 NK_WINDOW_BORDER | NK_WINDOW_MOVABLE | NK_WINDOW_SCALABLE |
                     NK_WINDOW_MINIMIZABLE | NK_WINDOW_TITLE)) {
        nk_layout_row_dynamic(nuklear_ctx, 30, 1);
        nk_label(nuklear_ctx, pars->filename.c_str(), NK_TEXT_LEFT);
        nk_layout_row_dynamic(nuklear_ctx, 30, 2);
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
        pars->srgb = nk_check_label(nuklear_ctx, "srgb output", pars->srgb);
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
    draw_widgets(window);
    glfwSwapBuffers(window);
}

void run_ui(yshade_app::params* pars) {
    // window
    auto window = yui::init_glfw(pars->width, pars->height, "yimview",
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

        // make ui
        draw_widgets(window);

        // swap buffers
        glfwSwapBuffers(window);

        // check for screenshot
        if (pars->no_ui) {
            save_screenshot(window, pars->imfilename);
            break;
        }

        // event hadling
        glfwWaitEvents();
    }

    yui::clear_nuklear((nk_context*)pars->widget_ctx, pars->legacy_gl);
    yui::clear_glfw(window);
}

int main(int argc, char* argv[]) {
    // params
    auto params = new yshade_app::params();
    auto parser = ycmd::make_parser(argc, argv, "trace meshes");
    yshade_app::init_params(params, parser);
    ycmd::check_parser(parser);

    // run ui
    run_ui(params);

    // done
    return EXIT_SUCCESS;
}
