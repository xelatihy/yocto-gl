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

#include "yapp.h"
#include "yui.h"

#include "../yocto/yocto_cmd.h"
#include "../yocto/yocto_math.h"
#include "../yocto/yocto_shape.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../yocto/stb_image_write.h"

void save_screenshot(const yui::info& info, const std::string& imfilename) {
    if (ycmd::get_extension(imfilename) != ".png") {
        printf("supports only png screenshots");
        return;
    }

    auto w = info.framebuffer_size[0];
    auto h = info.framebuffer_size[1];
    auto pixels = std::vector<unsigned char>(w * h * 4);
    glReadBuffer(GL_FRONT);
    glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, pixels.data());
    std::vector<unsigned char> line(w * 4);
    for (int j = 0; j < h / 2; j++) {
        memcpy(line.data(), pixels.data() + j * w * 4, w * 4);
        memcpy(pixels.data() + j * w * 4, pixels.data() + (h - 1 - j) * w * 4,
               w * 4);
        memcpy(pixels.data() + (h - 1 - j) * w * 4, line.data(), w * 4);
    }
    stbi_write_png(imfilename.c_str(), w, h, 4, pixels.data(), w * 4);
}

int main(int argc, char* argv[]) {
    // command line
    auto parser = ycmd::make_parser(argc, argv, "view meshes");
    auto exposure =
        ycmd::parse_opt<float>(parser, "--exposure", "-e", "image exposure", 0);
    auto gamma =
        ycmd::parse_opt<float>(parser, "--gamma", "-g", "image gamma", 2.2);
    auto amb =
        ycmd::parse_opt<float>(parser, "--ambient", "", "ambient factor", 0);
    auto camera_lights = ycmd::parse_flag(parser, "--camera_lights", "-c",
                                          "enable camera lights", false);
    auto camera = ycmd::parse_opt<int>(parser, "--camera", "-C", "camera", 0);
    auto no_ui = ycmd::parse_flag(parser, "--no-ui", "", "runs offline", false);
    auto aspect = ycmd::parse_opt<float>(parser, "--aspect", "-a",
                                         "image aspect", 16.0f / 9.0f);
    auto res = ycmd::parse_opt<int>(parser, "--resolution", "-r",
                                    "image resolution", 720);
    auto imfilename = ycmd::parse_opt<std::string>(parser, "--output", "-o",
                                                   "image filename", "out.png");
    auto filename = ycmd::parse_arg<std::string>(parser, "scene",
                                                 "scene filename", "", true);
    ycmd::check_parser(parser);

    // load scene
    auto scene = yapp::load_scene(filename);
    scene.cameras[camera].aspect = aspect;

    // view variables
    auto cur_background = 0;
    const std::array<ym::vec4f, 4> backgrounds = {{{0.0f, 0.0f, 0.0f, 0.0f},
                                                   {0.18f, 0.18f, 0.18f, 0.0f},
                                                   {0.5f, 0.5f, 0.5f, 0.0f},
                                                   {1.0f, 1.0f, 1.0f, 0.0f}}};
    auto wireframe = false;
    auto edges = false;

    // shading
    auto shade_prog = 0;
    std::vector<int> shade_txt;

    // prepare ui context
    auto context = yui::context();

    // init callback
    context.init += std::function<void(const yui::info& info)>(
        [&](const yui::info& info) {  // load textures
            yapp::init_shade(scene, shade_prog, shade_txt);
        });

    // window size callback
    context.window_size +=
        std::function<void(const yui::info& info)>([&](const yui::info& info) {
            auto& cam = scene.cameras[camera];
            cam.aspect = (float)info.win_size[0] / (float)info.win_size[1];
        });

    // window refresh callback
    context.window_refresh +=
        std::function<void(const yui::info& info)>([&](const yui::info& info) {
            // draw
            yapp::shade(scene, camera, shade_prog, shade_txt,
                        backgrounds[cur_background], exposure, gamma, wireframe,
                        edges, camera_lights);
        });

    // check continue callback
    context.update +=
        std::function<int(const yui::info& info)>([&](const yui::info& info) {
            // check for screenshot
            if (no_ui) {
                save_screenshot(info, imfilename);
                return -1;
            }

            // continue as usual
            return 0;
        });

    // text callback
    context.text += std::function<void(const yui::info& info, unsigned int)>(
        [&](const yui::info& info, unsigned int key) {
            switch (key) {
                case '[': exposure -= 1; break;
                case ']': exposure += 1; break;
                case '{': gamma -= 0.1f; break;
                case '}': gamma += 0.1f; break;
                case 'w': wireframe = !wireframe; break;
                case 'e': edges = !edges; break;
                case 'b':
                    cur_background = (cur_background + 1) % backgrounds.size();
                    break;
                case 's': save_screenshot(info, imfilename); break;
                case 'c': camera_lights = !camera_lights; break;
                case 'C': camera = (camera + 1) % scene.cameras.size(); break;
                case 't': {
                    for (auto& shape : scene.shapes) {
                        yshape::tesselate_stdshape(
                            shape.lines, shape.triangles, shape.pos, shape.norm,
                            shape.texcoord, shape.color, shape.radius);
                        int max_id = -1, min_id = 10000000;
                        for (auto& t : shape.triangles) {
                            for (auto i = 0; i < 3; i++)
                                max_id = ym::max(max_id, t[i]);
                            for (auto i = 0; i < 3; i++)
                                min_id = ym::min(min_id, t[i]);
                        }
                    }
                } break;
                default: printf("unsupported key\n"); break;
            }
        });

    // mouse position callback
    context.mouse_pos +=
        std::function<void(const yui::info& info)>([&](const yui::info& info) {
            if (info.mouse_button) {
                auto dolly = 0.0f;
                auto pan = ym::zero2f;
                auto rotate = ym::zero2f;
                switch (info.mouse_button) {
                    case 1:
                        rotate = (info.mouse_pos - info.mouse_last) / 100;
                        break;
                    case 2:
                        dolly =
                            (info.mouse_pos[0] - info.mouse_last[0]) / 100.0f;
                        break;
                    case 3:
                        pan = (info.mouse_pos - info.mouse_last) / 100;
                        break;
                    default: break;
                }

                auto& cam = scene.cameras[camera];
                ym::turntable(cam.frame, cam.focus, rotate, dolly, pan);
            }
        });

    // run ui
    yui::ui_loop(context, (int)std::round(aspect * res), res, "yview");

    // done
    return EXIT_SUCCESS;
}
