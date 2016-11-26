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
    std::string filename;
    std::vector<float> pixelf;
    std::vector<unsigned char> pixelb;
    int w, h, nc;
    int tex_glid;
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
    static const std::string labels[] = {"r", "g", "b", "a"};

    auto ij = ym::vec2i(round((pos - offset) / zoom));

    auto buf = std::string();
    buf += "img: " + std::to_string(img.w) + " x " + std::to_string(img.h) +
           " @ " + std::to_string(img.nc) + "\n";
    buf +=
        "mouse: " + std::to_string(ij[0]) + " " + std::to_string(ij[1]) + "\n";
    if (ij[0] >= 0 && ij[0] < img.w && ij[1] >= 0 && ij[1] < img.h) {
        auto cf = ym::vec4f(0, 0, 0, 1);
        auto cb = ym::vec4b(0, 0, 0, 255);
        if (img.pixelf.empty()) {
            for (auto i = 0; i < img.nc; i++) {
                cf[i] = img.pixelf[(ij[1] * img.w + ij[0]) * img.nc + i];
                cb[i] = ym::clamp(int(cf[i] * 256), 0, 255);
            }
        } else if (img.pixelb.empty()) {
            for (auto i = 0; i < img.nc; i++) {
                cb[i] = img.pixelb[(ij[1] * img.w + ij[0]) * img.nc + i];
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
}

std::vector<view_img> load_images(
    const std::vector<std::string>& img_filenames) {
    auto imgs = std::vector<view_img>();
    for (auto filename : img_filenames) {
        imgs.push_back(view_img());
        auto& img = imgs.back();
        img.filename = filename;
        auto ext = ycmd::get_extension(filename);
        if (ext == ".hdr") {
            auto pixels =
                stbi_loadf(filename.c_str(), &img.w, &img.h, &img.nc, 0);
            img.pixelf =
                std::vector<float>(pixels, pixels + img.w * img.h * img.nc);
            free(pixels);
        } else {
            auto pixels =
                stbi_load(filename.c_str(), &img.w, &img.h, &img.nc, 0);
            img.pixelb = std::vector<unsigned char>(
                pixels, pixels + img.w * img.h * img.nc);
            free(pixels);
        }
        if (img.pixelf.empty() && img.pixelb.empty()) {
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
    auto filenames = ycmd::parse_arga<std::string>(parser, "image",
                                                   "image filename", {}, true);
    ycmd::check_parser(parser);

    // loading images
    auto imgs = load_images(filenames);

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
    context.init += std::function<void(const yui::info& info)>(
        [&](const yui::info& info) {  // load textures
            for (auto& img : imgs) {
                if (!img.pixelf.empty()) {
                    img.tex_glid = yglu::make_texture(
                        img.w, img.h, img.nc, img.pixelf.data(), false, true);
                } else if (!img.pixelb.empty()) {
                    img.tex_glid = yglu::make_texture(
                        img.w, img.h, img.nc, img.pixelb.data(), false, true);
                } else
                    assert(false);
            }
        });

    // window refresh callback
    context.window_refresh +=
        std::function<void(const yui::info& info)>([&](const yui::info& info) {
            auto& img = imgs[cur_img];

            // begin frame
            auto background = backgrounds[cur_background];
            glClearColor(background[0], background[1], background[2],
                         background[3]);
            glDisable(GL_DEPTH_TEST);
            glClear(GL_COLOR_BUFFER_BIT);

            // draw image
            yglu::shade_image(img.tex_glid, img.w, img.h, img.w, img.h,
                              offset[0], offset[1], zoom, exposure, gamma);
        });

    // text callback
    context.text += std::function<void(const yui::info& info, unsigned int)>(
        [&](const yui::info& info, unsigned int key) {
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
        });

    // mouse position callback
    context.mouse_pos +=
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
        });

    // run ui
    yui::ui_loop(context, imgs[0].w, imgs[0].h, "yview");

    // done
    return EXIT_SUCCESS;
}
