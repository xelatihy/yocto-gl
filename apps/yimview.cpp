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

#include "yui.h"

#include "../yocto/yocto_cmd.h"
#include "../yocto/yocto_glu.h"
#include "../yocto/yocto_math.h"

#define STB_IMAGE_IMPLEMENTATION
#include "../yocto/stb_image.h"

struct app_img {
    // image path
    std::string filename;

    // original image data size
    int width = 0;
    int height = 0;
    int ncomp = 0;

    // pixel data in RGBA format
    ym::image<ym::vec4f> hdr;
    ym::image<ym::vec4b> ldr;

    // opengl texture
    yglu::uint tex_glid = 0;

    // hdr controls
    float hdr_exposure = 0;
    float hdr_gamma = 2.2f;

    // check hdr
    bool is_hdr() const { return !hdr.empty(); }
};

// images
std::vector<std::string> filenames;
std::vector<app_img> imgs;

// hdr controls
float hdr_exposure = 0;
float hdr_gamma = 2.2;

// opengl type
bool legacy_gl = false;

// view parameters
int cur_img = 0;
int cur_background = 0;
float zoom = 1;
ym::vec2f offset = ym::vec2f();
bool hud_print = false;
const std::array<ym::vec4f, 4> backgrounds = {{{0.0f, 0.0f, 0.0f, 0.0f},
                                               {0.18f, 0.18f, 0.18f, 0.0f},
                                               {0.5f, 0.5f, 0.5f, 0.0f},
                                               {1.0f, 1.0f, 1.0f, 0.0f}}};

// glfw
GLFWwindow* window = nullptr;

// nuklear
nk_context* nuklear_ctx = nullptr;
int hud_width = 256;

std::vector<app_img> load_images(const std::vector<std::string>& img_filenames,
                                 float exposure, float gamma) {
    auto imgs = std::vector<app_img>();
    for (auto filename : img_filenames) {
        imgs.push_back(app_img());
        auto& img = imgs.back();
        img.filename = filename;
        auto ext = ycmd::get_extension(filename);
        if (ext == ".hdr") {
            auto pixels = stbi_loadf(filename.c_str(), &img.width, &img.height,
                                     &img.ncomp, 0);
            img.hdr =
                ym::make_image4(img.width, img.height, img.ncomp, pixels, 1.0f);
            img.ldr.resize(img.hdr.size());
            ym::exposure_gamma(img.hdr, img.ldr, exposure, gamma);
            img.hdr_exposure = exposure;
            img.hdr_gamma = gamma;
            free(pixels);
        } else {
            auto pixels = stbi_load(filename.c_str(), &img.width, &img.height,
                                    &img.ncomp, 0);
            img.ldr = ym::make_image4(img.width, img.height, img.ncomp, pixels,
                                      (unsigned char)255);
            free(pixels);
        }
        if (img.hdr.empty() && img.ldr.empty()) {
            printf("cannot load image %s\n", img.filename.c_str());
            exit(1);
        }
        img.tex_glid = 0;
    }
    return imgs;
}

void text_callback(GLFWwindow* window, unsigned int key) {
    nk_glfw3_gl3_char_callback(window, key);
    if (nk_item_is_any_active(nuklear_ctx)) return;
    switch (key) {
        case ' ':
        case '.': cur_img = (cur_img + 1) % imgs.size(); break;
        case ',':
            cur_img = (cur_img - 1 + (int)imgs.size()) % imgs.size();
            break;
        case '-':
        case '_': zoom /= 2; break;
        case '+':
        case '=': zoom *= 2; break;
        case '[': hdr_exposure -= 1; break;
        case ']': hdr_exposure += 1; break;
        case '{': hdr_gamma -= 0.1f; break;
        case '}': hdr_gamma += 0.1f; break;
        case '1':
            hdr_exposure = 0;
            hdr_gamma = 1;
            break;
        case '2':
            hdr_exposure = 0;
            hdr_gamma = 2.2f;
            break;
        case 'z': zoom = 1; break;
        case 'b':
            cur_background = (cur_background + 1) % backgrounds.size();
            break;
        case 'h':
            // TODO: hud
            break;
        default: printf("unsupported key\n"); break;
    }
}

void draw_image() {
    auto framebuffer_size = yui::framebuffer_size(window);
    glViewport(0, 0, framebuffer_size[0], framebuffer_size[1]);

    auto& img = imgs[cur_img];

    // begin frame
    auto background = backgrounds[cur_background];
    glClearColor(background[0], background[1], background[2], background[3]);
    glDisable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT);

    // draw image
    auto window_size = yui::window_size(window);
    if (legacy_gl) {
        yglu::legacy::draw_image(img.tex_glid, img.width, img.height,
                                 window_size[0], window_size[1], offset[0],
                                 offset[1], zoom);
    } else {
        yglu::modern::shade_image(img.tex_glid, img.width, img.height,
                                  window_size[0], window_size[1], offset[0],
                                  offset[1], zoom);
    }
}

void draw_widgets() {
    auto& img = imgs[cur_img];
    auto mouse_pos = yui::mouse_pos(window);
    auto window_size = yui::window_size(window);
    if (legacy_gl) {
        nk_glfw3_gl2_new_frame();
    } else {
        nk_glfw3_gl3_new_frame();
    }
    if (nk_begin(nuklear_ctx, "yimview", nk_rect(window_size[0] - hud_width, 0,
                                                 hud_width, window_size[1]),
                 NK_WINDOW_BORDER | NK_WINDOW_MOVABLE | NK_WINDOW_SCALABLE |
                     NK_WINDOW_MINIMIZABLE | NK_WINDOW_TITLE)) {
        nk_layout_row_dynamic(nuklear_ctx, 30, 1);
        nk_label(nuklear_ctx, img.filename.c_str(), NK_TEXT_LEFT);
        nk_layout_row_dynamic(nuklear_ctx, 30, 3);
        nk_value_int(nuklear_ctx, "w", img.width);
        nk_value_int(nuklear_ctx, "h", img.height);
        nk_value_int(nuklear_ctx, "c", img.ncomp);
        auto xy = (mouse_pos - offset) / zoom;
        auto ij = ym::vec2i(round(xy[0]), round(xy[1]));
        auto inside =
            ij[0] >= 0 && ij[1] >= 0 && ij[0] < img.width && ij[1] < img.height;
        nk_layout_row_dynamic(nuklear_ctx, 30, 4);
        nk_value_int(nuklear_ctx, "r", (inside) ? img.ldr[ij][0] : 0);
        nk_value_int(nuklear_ctx, "g", (inside) ? img.ldr[ij][1] : 0);
        nk_value_int(nuklear_ctx, "b", (inside) ? img.ldr[ij][2] : 0);
        nk_value_int(nuklear_ctx, "a", (inside) ? img.ldr[ij][3] : 0);
        if (img.is_hdr()) {
            nk_layout_row_dynamic(nuklear_ctx, 30, 2);
            nk_value_float(nuklear_ctx, "r", (inside) ? img.hdr[ij][0] : 0);
            nk_value_float(nuklear_ctx, "g", (inside) ? img.hdr[ij][1] : 0);
            nk_value_float(nuklear_ctx, "b", (inside) ? img.hdr[ij][2] : 0);
            nk_value_float(nuklear_ctx, "a", (inside) ? img.hdr[ij][3] : 0);
            nk_layout_row_dynamic(nuklear_ctx, 30, 1);
            nk_property_float(nuklear_ctx, "exposure", -20, &hdr_exposure, 20,
                              1, 1);
            nk_property_float(nuklear_ctx, "gamma", 0.1, &hdr_gamma, 5, 0.1,
                              0.1);
        }
    }
    nk_end(nuklear_ctx);

    if (legacy_gl) {
        nk_glfw3_gl2_render(NK_ANTI_ALIASING_ON, 512 * 1024, 128 * 1024);
    } else {
        nk_glfw3_gl3_render(NK_ANTI_ALIASING_ON, 512 * 1024, 128 * 1024);
    }
}

void window_refresh_callback(GLFWwindow* window) {
    draw_image();
    draw_widgets();
    glfwSwapBuffers(window);
}

void run_ui() {
    // window
    window = yui::init_glfw({imgs[0].width + hud_width, imgs[0].height},
                            "yimview", legacy_gl, nullptr, text_callback);

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

    nuklear_ctx = yui::init_nuklear(window, legacy_gl);

    // load textures
    for (auto& img : imgs) {
        if (legacy_gl) {
            img.tex_glid = yglu::legacy::make_texture(
                img.width, img.height, 4, (unsigned char*)img.ldr.data(), false,
                false);
        } else {
            img.tex_glid = yglu::modern::make_texture(
                img.width, img.height, 4, (unsigned char*)img.ldr.data(), false,
                false, false);
        }
    }

    while (!glfwWindowShouldClose(window)) {
        glfwGetWindowSize(window, &window_size[0], &window_size[1]);
        glfwGetFramebufferSize(window, &framebuffer_size[0],
                               &framebuffer_size[1]);

        mouse_last = mouse_pos;
        mouse_pos = yui::mouse_pos(window);
        mouse_button = yui::mouse_button(window);

        auto& img = imgs[cur_img];
        glfwSetWindowTitle(
            window,
            ("yimview | " + img.filename + " | " + std::to_string(img.width) +
             "x" + std::to_string(img.height) + "@" + std::to_string(img.ncomp))
                .c_str());

        // handle mouse
        if (mouse_button && mouse_pos != mouse_last &&
            !nk_item_is_any_active(nuklear_ctx)) {
            switch (mouse_button) {
                case 1: offset += mouse_pos - mouse_last; break;
                case 2:
                    zoom *= powf(2, (mouse_pos[0] - mouse_last[0]) * 0.001f);
                    break;
                default: break;
            }
        }

        // refresh hdr
        if (img.is_hdr() &&
            (hdr_exposure != img.hdr_exposure || hdr_gamma != img.hdr_gamma)) {
            ym::exposure_gamma(img.hdr, img.ldr, hdr_exposure, hdr_gamma);
            img.hdr_exposure = hdr_exposure;
            img.hdr_gamma = hdr_gamma;
            if (legacy_gl) {
                yglu::legacy::update_texture(img.tex_glid, img.width,
                                             img.height, 4,
                                             img.ldr.data()->data(), false);
            } else {
                yglu::modern::update_texture(img.tex_glid, img.width,
                                             img.height, 4,
                                             img.ldr.data()->data(), false);
            }
        }

        // draw
        draw_image();

        // make ui
        draw_widgets();

        // swap buffers
        glfwSwapBuffers(window);

        // event hadling
        glfwWaitEvents();
    }

    yui::clear_nuklear(nuklear_ctx, legacy_gl);
    yui::clear_glfw(window);
}

int main(int argc, char* argv[]) {
    // command line params
    auto parser = ycmd::make_parser(argc, argv, "view images");
    hdr_exposure =
        ycmd::parse_opt<float>(parser, "--exposure", "-e", "image exposure", 0);
    hdr_gamma =
        ycmd::parse_opt<float>(parser, "--gamma", "-g", "image gamma", 2.2);
    legacy_gl = ycmd::parse_flag(parser, "--legacy_opengl", "-L",
                                 "uses legacy OpenGL", false);
    filenames = ycmd::parse_arga<std::string>(parser, "image", "image filename",
                                              {}, true);
    ycmd::check_parser(parser);

    // loading images
    imgs = load_images(filenames, hdr_exposure, hdr_gamma);

    // run ui
    run_ui();

    // done
    return EXIT_SUCCESS;
}
