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

struct view_img {
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
    int tex_glid = 0;

    // hdr controls
    float hdr_exposure = 0;
    float hdr_gamma = 2.2f;

    // check hdr
    bool is_hdr() const { return !hdr.empty(); }
};

// inspect string
#define strcatf(str, fmt, ...)                                                 \
    {                                                                          \
        char __buf[512];                                                       \
        sprintf(__buf, fmt, __VA_ARGS__);                                      \
        strcat(str, __buf);                                                    \
    }

std::string inspect_str(const view_img& img, const ym::vec2f& offset,
                        float zoom, const ym::vec2f& pos, float exposure,
                        float gamma) {
#if 0
    static const std::string labels[] = {"r", "g", "b", "a"};

    auto xy = (pos - offset) / zoom;
    auto ij = ym::vec2i(round(xy[0]), round(xy[1]));

    auto buf = std::string();
    buf += "img: " + std::to_string(img.width) + " x " +
                    std::to_string(img.height) +
           " @ " + std::to_string(img.ncomp) + "\n";
    buf +=
        "mouse: " + std::to_string(ij[0]) + " " + std::to_string(ij[1]) + "\n";
    if (ij[0] >= 0 && ij[0] < img.width && ij[1] >= 0 && ij[1] < img.height) {
        auto cf = ym::vec4f();
        auto cb = ym::vec4b();
        if (img.hdr.empty()) {
            for (auto i = 0; i < img.ncomp; i++) {
                cf[i] = img.hdr[(ij[1] * img.width + ij[0]) * img.ncomp + i];
                cb[i] = ym::clamp(int(cf[i] * 256), 0, 255);
            }
        } else if (img.ldr.empty()) {
            for (auto i = 0; i < img.ncomp; i++) {
                cb[i] = img.ldr[(ij[1] * img.width + ij[0]) * img.ncomp + i];
                cf[i] = std::pow(cb[i] / 255.0f, 1 / 2.2f);
            }
        } else
            assert(false);
        for (auto i = 0; i < 4; i++) {
            buf += labels[i] + ": " + std::to_string(cb[i]) + " " +
                   std::to_string(cf[i]) + "\n";
        }
    } else {
        for (auto i = 0; i < 4; i++) buf += labels[i] + ": -\n";
    }
    return buf;
#endif
    return {};
}

std::vector<view_img> load_images(const std::vector<std::string>& img_filenames,
                                  float exposure, float gamma) {
    auto imgs = std::vector<view_img>();
    for (auto filename : img_filenames) {
        imgs.push_back(view_img());
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

int main(int argc, char* argv[]) {
    // command line params
    auto parser = ycmd::make_parser(argc, argv, "view images");
    auto exposure =
        ycmd::parse_opt<float>(parser, "--exposure", "-e", "image exposure", 0);
    auto gamma =
        ycmd::parse_opt<float>(parser, "--gamma", "-g", "image gamma", 2.2);
    auto legacy_gl = ycmd::parse_flag(parser, "--legacy_opengl", "-L",
                                      "uses legacy OpenGL", false);
    auto filenames = ycmd::parse_arga<std::string>(parser, "image",
                                                   "image filename", {}, true);
    ycmd::check_parser(parser);

    // loading images
    auto imgs = load_images(filenames, exposure, gamma);

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

    // prepare ui context
    auto context = yui::context();

    // init callback
    context.init.push_back(std::function<void(const yui::info& info)>(
        [&](const yui::info& info) {  // load textures
            for (auto& img : imgs) {
                if (legacy_gl) {
                    img.tex_glid = yglu::legacy::make_texture(
                        img.width, img.height, 4,
                        (unsigned char*)img.ldr.data(), false, false);
                } else {
                    img.tex_glid = yglu::modern::make_texture(
                        img.width, img.height, 4,
                        (unsigned char*)img.ldr.data(), false, false, false);
                }
            }
        }));

    // window refresh callback
    context.window_refresh.push_back(
        std::function<void(const yui::info& info)>([&](const yui::info& info) {
            auto& img = imgs[cur_img];

            // refresh hdr
            if (img.is_hdr() &&
                (exposure != img.hdr_exposure || gamma != img.hdr_gamma)) {
                ym::exposure_gamma(img.hdr, img.ldr, exposure, gamma);
                img.hdr_exposure = exposure;
                img.hdr_gamma = gamma;
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

            // begin frame
            auto background = backgrounds[cur_background];
            glClearColor(background[0], background[1], background[2],
                         background[3]);
            glDisable(GL_DEPTH_TEST);
            glClear(GL_COLOR_BUFFER_BIT);

            // draw image
            if (legacy_gl) {
                yglu::legacy::draw_image(img.tex_glid, img.width, img.height,
                                         info.win_size[0], info.win_size[1],
                                         offset[0], offset[1], zoom);
            } else {
                yglu::modern::shade_image(img.tex_glid, img.width, img.height,
                                          info.win_size[0], info.win_size[1],
                                          offset[0], offset[1], zoom);
            }
        }));

    // text callback
    context.text.push_back(
        std::function<void(const yui::info& info, unsigned int)>([&](
            const yui::info& info, unsigned int key) {
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
                case '[': exposure -= 1; break;
                case ']': exposure += 1; break;
                case '{': gamma -= 0.1f; break;
                case '}': gamma += 0.1f; break;
                case '1':
                    exposure = 0;
                    gamma = 1;
                    break;
                case '2':
                    exposure = 0;
                    gamma = 2.2f;
                    break;
                case 'z': zoom = 1; break;
                case 'b':
                    cur_background = (cur_background + 1) % backgrounds.size();
                    break;
                case 'h':
                    // TODO: hud
                    break;
                case 'p': hud_print = !hud_print; break;
                default: printf("unsupported key\n"); break;
            }
        }));

    // mouse position callback
    context.mouse_pos.push_back(
        std::function<void(const yui::info& info)>([&](const yui::info& info) {
            switch (info.mouse_button) {
                case 1: offset += info.mouse_pos - info.mouse_last; break;
                case 2:
                    zoom *= powf(2, (info.mouse_pos[0] - info.mouse_last[0]) *
                                        0.001f);
                    break;
                default: break;
            }

            // draw inspector
            if (hud_print) {
                auto str = inspect_str(imgs[cur_img], offset, zoom,
                                       info.mouse_pos, exposure, gamma);
                printf("\033[2J\033[1;1H");
                printf("%s", str.c_str());
            }
        }));

    // run ui
    yui::ui_loop(context, imgs[0].width, imgs[0].height, "yimview", !legacy_gl);

    // done
    return EXIT_SUCCESS;
}
