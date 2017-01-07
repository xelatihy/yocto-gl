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

#include "yui.h"

#include "../yocto/yocto_cmd.h"
#include "../yocto/yocto_glu.h"
#include "../yocto/yocto_math.h"

#define STB_IMAGE_IMPLEMENTATION
#include "../yocto/stb_image.h"

namespace yimview_app {

struct img {
    // image path
    std::string filename;

    // original image data size
    int width = 0;
    int height = 0;
    int ncomp = 0;

    // pixel data
    std::vector<float> hdr;
    std::vector<unsigned char> ldr;

    // opengl texture
    yglu::uint tex_glid = 0;

    // hdr controls
    float exposure = 0;
    float gamma = 2.2f;
    bool srgb = true;

    // check hdr
    bool is_hdr() const { return !hdr.empty(); }
};

struct params {
    std::vector<std::string> filenames;
    std::vector<img*> imgs;

    float exposure = 0;
    float gamma = 1;
    bool srgb = true;

    bool legacy_gl = false;

    int cur_img = 0;
    int cur_background = 0;
    float zoom = 1;
    ym::vec2f offset = ym::vec2f();

    float background = 0;

    void* widget_ctx = nullptr;

    ~params() {
        for (auto img : imgs) delete img;
    }
};

std::vector<img*> load_images(const std::vector<std::string>& img_filenames,
                              float exposure, float gamma, bool srgb) {
    auto imgs = std::vector<img*>();
    for (auto filename : img_filenames) {
        imgs.push_back(new img());
        auto img = imgs.back();
        img->filename = filename;
        auto ext = ycmd::get_extension(filename);
        if (ext == ".hdr") {
            auto pixels = stbi_loadf(filename.c_str(), &img->width,
                                     &img->height, &img->ncomp, 0);
            img->hdr = std::vector<float>(
                pixels, pixels + img->width * img->height * img->ncomp);
            img->ldr.resize(img->hdr.size());
            ym::exposure_gamma(img->width, img->height, img->ncomp,
                               img->hdr.data(), img->ldr.data(), exposure,
                               gamma, srgb);
            img->exposure = exposure;
            img->gamma = gamma;
            free(pixels);
        } else {
            auto pixels = stbi_load(filename.c_str(), &img->width, &img->height,
                                    &img->ncomp, 0);
            img->ldr = std::vector<unsigned char>(
                pixels, pixels + img->width * img->height * img->ncomp);
            free(pixels);
        }
        if (img->hdr.empty() && img->ldr.empty()) {
            printf("cannot load image %s\n", img->filename.c_str());
            exit(1);
        }
        img->tex_glid = 0;
    }
    return imgs;
}

void init_params(params* pars, ycmd::parser* parser) {
    pars->exposure = ycmd::parse_opt<float>(parser, "--exposure", "-e",
                                            "hdr image exposure", 0);
    pars->gamma =
        ycmd::parse_opt<float>(parser, "--gamma", "-g", "hdr image gamma", 1);
    pars->srgb = ycmd::parse_opt<bool>(parser, "--srgb", "",
                                       "hdr image srgb output", true);
    pars->legacy_gl = ycmd::parse_flag(parser, "--legacy_opengl", "-L",
                                       "uses legacy OpenGL", false);
    auto filenames = ycmd::parse_arga<std::string>(parser, "image",
                                                   "image filename", {}, true);

    // loading images
    pars->imgs =
        load_images(filenames, pars->exposure, pars->gamma, pars->srgb);
}
}  // namespace

const int hud_width = 256;

void text_callback(GLFWwindow* window, unsigned int key) {
    auto pars = (yimview_app::params*)glfwGetWindowUserPointer(window);
    auto nuklear_ctx = (nk_context*)pars->widget_ctx;
    nk_glfw3_gl3_char_callback(window, key);
    if (nk_item_is_any_active(nuklear_ctx)) return;
    switch (key) {
        case ' ':
        case '.':
            pars->cur_img = (pars->cur_img + 1) % pars->imgs.size();
            break;
        case ',':
            pars->cur_img = (pars->cur_img - 1 + (int)pars->imgs.size()) %
                            pars->imgs.size();
            break;
        case '-':
        case '_': pars->zoom /= 2; break;
        case '+':
        case '=': pars->zoom *= 2; break;
        case '[': pars->exposure -= 1; break;
        case ']': pars->exposure += 1; break;
        case '{': pars->gamma -= 0.1f; break;
        case '}': pars->gamma += 0.1f; break;
        case '\\': pars->srgb = pars->srgb; break;
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
        case 'z': pars->zoom = 1; break;
        case 'h':
            // TODO: hud
            break;
        default: printf("unsupported key\n"); break;
    }
}

void draw_image(GLFWwindow* window) {
    auto pars = (yimview_app::params*)glfwGetWindowUserPointer(window);
    auto framebuffer_size = yui::framebuffer_size(window);
    glViewport(0, 0, framebuffer_size[0], framebuffer_size[1]);

    auto img = pars->imgs[pars->cur_img];

    // begin frame
    glClearColor(pars->background, pars->background, pars->background, 0);
    glDisable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT);

    // draw image
    auto window_size = yui::window_size(window);
    if (pars->legacy_gl) {
        yglu::legacy::draw_image(img->tex_glid, img->width, img->height,
                                 window_size[0], window_size[1],
                                 pars->offset[0], pars->offset[1], pars->zoom);
    } else {
        yglu::modern::shade_image(img->tex_glid, img->width, img->height,
                                  window_size[0], window_size[1],
                                  pars->offset[0], pars->offset[1], pars->zoom);
    }
}

template <typename T>
ym::vec<T, 4> lookup_image(int w, int h, int nc, const T* pixels, int x, int y,
                           T one) {
    if (x < 0 || y < 0 || x > w - 1 || y > h - 1) return {0, 0, 0, 0};
    auto v = ym::vec<T, 4>{0, 0, 0, 0};
    auto vv = pixels + ((w * y) + x) * nc;
    switch (nc) {
        case 1: v = {vv[0], 0, 0, one}; break;
        case 2: v = {vv[0], vv[1], 0, one}; break;
        case 3: v = {vv[0], vv[1], vv[2], one}; break;
        case 4: v = {vv[0], vv[1], vv[2], vv[3]}; break;
        default: assert(false);
    }
    return v;
}

void draw_widgets(GLFWwindow* window) {
    auto pars = (yimview_app::params*)glfwGetWindowUserPointer(window);
    auto nuklear_ctx = (nk_context*)pars->widget_ctx;
    auto& img = pars->imgs[pars->cur_img];
    auto mouse_pos = yui::mouse_pos(window);
    auto window_size = yui::window_size(window);
    if (pars->legacy_gl) {
        nk_glfw3_gl2_new_frame();
    } else {
        nk_glfw3_gl3_new_frame();
    }
    if (nk_begin(nuklear_ctx, "yimview", nk_rect(window_size[0] - hud_width, 0,
                                                 hud_width, window_size[1]),
                 NK_WINDOW_BORDER | NK_WINDOW_MOVABLE | NK_WINDOW_SCALABLE |
                     NK_WINDOW_MINIMIZABLE | NK_WINDOW_TITLE)) {
        nk_layout_row_dynamic(nuklear_ctx, 30, 1);
        nk_label(nuklear_ctx, img->filename.c_str(), NK_TEXT_LEFT);
        nk_layout_row_dynamic(nuklear_ctx, 30, 3);
        nk_value_int(nuklear_ctx, "w", img->width);
        nk_value_int(nuklear_ctx, "h", img->height);
        nk_value_int(nuklear_ctx, "c", img->ncomp);
        auto xy = (mouse_pos - pars->offset) / pars->zoom;
        auto ij = ym::vec2i{(int)round(xy[0]), (int)round(xy[1])};
        auto inside = ij[0] >= 0 && ij[1] >= 0 && ij[0] < img->width &&
                      ij[1] < img->height;
        nk_layout_row_dynamic(nuklear_ctx, 30, 4);
        auto ldrp =
            lookup_image(img->width, img->height, img->ncomp, img->ldr.data(),
                         ij[0], ij[1], (unsigned char)255);
        nk_value_int(nuklear_ctx, "r", (inside) ? ldrp[0] : 0);
        nk_value_int(nuklear_ctx, "g", (inside) ? ldrp[1] : 0);
        nk_value_int(nuklear_ctx, "b", (inside) ? ldrp[2] : 0);
        nk_value_int(nuklear_ctx, "a", (inside) ? ldrp[3] : 0);
        if (img->is_hdr()) {
            auto hdrp = lookup_image(img->width, img->height, img->ncomp,
                                     img->hdr.data(), ij[0], ij[1], 1.0f);
            nk_layout_row_dynamic(nuklear_ctx, 30, 2);
            nk_value_float(nuklear_ctx, "r", (inside) ? hdrp[0] : 0);
            nk_value_float(nuklear_ctx, "g", (inside) ? hdrp[1] : 0);
            nk_value_float(nuklear_ctx, "b", (inside) ? hdrp[2] : 0);
            nk_value_float(nuklear_ctx, "a", (inside) ? hdrp[3] : 0);
            nk_layout_row_dynamic(nuklear_ctx, 30, 1);
            nk_property_float(nuklear_ctx, "hdr exposure", -20, &pars->exposure,
                              20, 1, 1);
            nk_property_float(nuklear_ctx, "hdr gamma", 0.1, &pars->gamma, 5,
                              0.1, 0.1);
            pars->srgb =
                nk_check_label(nuklear_ctx, "hdr srgb output", pars->srgb);
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
    draw_image(window);
    draw_widgets(window);
    glfwSwapBuffers(window);
}

void run_ui(yimview_app::params* pars) {
    // window
    auto window =
        yui::init_glfw(pars->imgs[0]->width + hud_width, pars->imgs[0]->height,
                       "yimview", pars->legacy_gl, pars, text_callback);

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

    // load textures
    for (auto& img : pars->imgs) {
        if (pars->legacy_gl) {
            img->tex_glid = yglu::legacy::make_texture(
                img->width, img->height, img->ncomp,
                (unsigned char*)img->ldr.data(), false, false);
        } else {
            img->tex_glid = yglu::modern::make_texture(
                img->width, img->height, img->ncomp,
                (unsigned char*)img->ldr.data(), false, false, false);
        }
    }

    while (!glfwWindowShouldClose(window)) {
        glfwGetWindowSize(window, &window_size[0], &window_size[1]);
        glfwGetFramebufferSize(window, &framebuffer_size[0],
                               &framebuffer_size[1]);

        mouse_last = mouse_pos;
        mouse_pos = yui::mouse_pos(window);
        mouse_button = yui::mouse_button(window);

        auto& img = pars->imgs[pars->cur_img];
        glfwSetWindowTitle(window, ("yimview | " + img->filename + " | " +
                                    std::to_string(img->width) + "x" +
                                    std::to_string(img->height) + "@" +
                                    std::to_string(img->ncomp))
                                       .c_str());

        // handle mouse
        if (mouse_button && mouse_pos != mouse_last &&
            !nk_item_is_any_active((nk_context*)pars->widget_ctx)) {
            switch (mouse_button) {
                case 1: pars->offset += mouse_pos - mouse_last; break;
                case 2:
                    pars->zoom *=
                        powf(2, (mouse_pos[0] - mouse_last[0]) * 0.001f);
                    break;
                default: break;
            }
        }

        // refresh hdr
        if (img->is_hdr() &&
            (pars->exposure != img->exposure || pars->gamma != img->gamma ||
             pars->srgb != img->srgb)) {
            ym::exposure_gamma(img->width, img->height, img->ncomp,
                               img->hdr.data(), img->ldr.data(), pars->exposure,
                               pars->gamma, pars->srgb);
            img->exposure = pars->exposure;
            img->gamma = pars->gamma;
            img->srgb = pars->srgb;
            if (pars->legacy_gl) {
                yglu::legacy::update_texture(img->tex_glid, img->width,
                                             img->height, img->ncomp,
                                             img->ldr.data(), false);
            } else {
                yglu::modern::update_texture(img->tex_glid, img->width,
                                             img->height, img->ncomp,
                                             img->ldr.data(), false);
            }
        }

        // draw
        draw_image(window);
        draw_widgets(window);

        // swap buffers
        glfwSwapBuffers(window);

        // event hadling
        glfwWaitEvents();
    }

    yui::clear_nuklear((nk_context*)pars->widget_ctx, pars->legacy_gl);
    yui::clear_glfw(window);
}

int main(int argc, char* argv[]) {
    // command line params
    auto pars = new yimview_app::params();
    auto parser = ycmd::make_parser(argc, argv, "view images");
    yimview_app::init_params(pars, parser);
    ycmd::check_parser(parser);

    // run ui
    run_ui(pars);

    // done
    delete pars;
    return EXIT_SUCCESS;
}
