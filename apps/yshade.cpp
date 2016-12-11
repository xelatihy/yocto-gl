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

// scene
std::string filename;
std::string imfilename;
yapp::scene scene;

// lighting
float hdr_exposure = 0;
float hdr_gamma = 2.2;
float amb = 0;
bool camera_lights = false;

// camera
int camera = 0;
float aspect = 16.0f / 9.0f;
int res = 720;

// gl
bool no_ui = false;
bool legacy_gl = false;

// view variables
int cur_background = 0;
const std::array<ym::vec4f, 4> backgrounds = {{{0.0f, 0.0f, 0.0f, 0.0f},
                                               {0.18f, 0.18f, 0.18f, 0.0f},
                                               {0.5f, 0.5f, 0.5f, 0.0f},
                                               {1.0f, 1.0f, 1.0f, 0.0f}}};
bool wireframe = false;
bool edges = false;

// shading
yglu::uint shade_prog = 0;
yglu::uint shade_vao = 0;
std::vector<yglu::uint> shade_txt;
std::vector<std::array<yglu::uint, 7>> shade_vbo;

// glfw
GLFWwindow* window = nullptr;

// nuklear
nk_context* nuklear_ctx = nullptr;
int hud_width = 256;

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

void draw_scene() {
    auto window_size = yui::window_size(window);
    scene.cameras[camera].aspect =
        (float)window_size[0] / (float)window_size[1];
    if (legacy_gl) {
        yapp::draw(scene, camera, shade_txt, backgrounds[cur_background],
                   hdr_exposure, hdr_gamma, wireframe, edges, camera_lights,
                   {amb, amb, amb});
    } else {
        yapp::shade(scene, camera, shade_prog, shade_vao, shade_txt, shade_vbo,
                    backgrounds[cur_background], hdr_exposure, hdr_gamma,
                    wireframe, edges, camera_lights, {amb, amb, amb});
    }
}

void text_callback(GLFWwindow* window, unsigned int key) {
    nk_glfw3_gl3_char_callback(window, key);
    if (nk_item_is_any_active(nuklear_ctx)) return;
    switch (key) {
        case '[': hdr_exposure -= 1; break;
        case ']': hdr_exposure += 1; break;
        case '{': hdr_gamma -= 0.1f; break;
        case '}': hdr_gamma += 0.1f; break;
        case 'w': wireframe = !wireframe; break;
        case 'e': edges = !edges; break;
        case 'b':
            cur_background = (cur_background + 1) % backgrounds.size();
            break;
        case 's': save_screenshot(window, imfilename); break;
        case 'c': camera_lights = !camera_lights; break;
        case 'C': camera = (camera + 1) % scene.cameras.size(); break;
        case 't': {
            for (auto& shape : scene.shapes) {
                yshape::tesselate_stdshape(
                    shape.lines, shape.triangles, shape.pos, shape.norm,
                    shape.texcoord, shape.color, shape.radius);
            }
        } break;
        default: printf("unsupported key\n"); break;
    }
}

void draw_widgets() {
    auto window_size = yui::window_size(window);
    if (legacy_gl) {
        nk_glfw3_gl2_new_frame();
    } else {
        nk_glfw3_gl3_new_frame();
    }
    if (nk_begin(nuklear_ctx, "yshade", nk_rect(window_size[0] - hud_width, 0,
                                                hud_width, window_size[1]),
                 NK_WINDOW_BORDER | NK_WINDOW_MOVABLE | NK_WINDOW_SCALABLE |
                     NK_WINDOW_MINIMIZABLE | NK_WINDOW_TITLE)) {
        nk_layout_row_dynamic(nuklear_ctx, 30, 1);
        nk_label(nuklear_ctx, filename.c_str(), NK_TEXT_LEFT);
        nk_layout_row_dynamic(nuklear_ctx, 30, 2);
        nk_property_int(nuklear_ctx, "camera", 0, &camera,
                        (int)scene.cameras.size() - 1, 1, 1);
        camera_lights = nk_check_label(nuklear_ctx, "eyelight", camera_lights);
        wireframe = nk_check_label(nuklear_ctx, "wireframe", wireframe);
        edges = nk_check_label(nuklear_ctx, "edges", edges);
        nk_layout_row_dynamic(nuklear_ctx, 30, 1);
        nk_property_float(nuklear_ctx, "exposure", -20, &hdr_exposure, 20, 1,
                          1);
        nk_property_float(nuklear_ctx, "gamma", 0.1, &hdr_gamma, 5, 0.1, 0.1);
        if (nk_button_label(nuklear_ctx, "tesselate")) {
            for (auto& shape : scene.shapes) {
                yshape::tesselate_stdshape(
                    shape.lines, shape.triangles, shape.pos, shape.norm,
                    shape.texcoord, shape.color, shape.radius);
            }
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
    draw_scene();
    draw_widgets();
    glfwSwapBuffers(window);
}

void run_ui() {
    // window
    window = yui::init_glfw({(int)(aspect * res), res}, "yimview", legacy_gl,
                            nullptr, text_callback);

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
    if (legacy_gl) {
        yapp::init_draw(scene, shade_txt);
    } else {
        yapp::init_shade(scene, shade_prog, shade_vao, shade_txt, shade_vbo);
    }

    nuklear_ctx = yui::init_nuklear(window, legacy_gl);

    while (!glfwWindowShouldClose(window)) {
        glfwGetWindowSize(window, &window_size[0], &window_size[1]);
        glfwGetFramebufferSize(window, &framebuffer_size[0],
                               &framebuffer_size[1]);

        mouse_last = mouse_pos;
        mouse_pos = yui::mouse_pos(window);
        mouse_button = yui::mouse_button(window);

        glfwSetWindowTitle(window, ("yshade | " + filename).c_str());

        // handle mouse
        if (mouse_button && mouse_pos != mouse_last &&
            !nk_item_is_any_active(nuklear_ctx)) {
            auto dolly = 0.0f;
            auto pan = ym::zero2f;
            auto rotate = ym::zero2f;
            switch (mouse_button) {
                case 1: rotate = (mouse_pos - mouse_last) / 100; break;
                case 2: dolly = (mouse_pos[0] - mouse_last[0]) / 100.0f; break;
                case 3: pan = (mouse_pos - mouse_last) / 100; break;
                default: break;
            }

            auto& cam = scene.cameras[camera];
            ym::turntable(cam.frame, cam.focus, rotate, dolly, pan);
        }

        // draw
        draw_scene();

        // make ui
        draw_widgets();

        // swap buffers
        glfwSwapBuffers(window);

        // check for screenshot
        if (no_ui) {
            save_screenshot(window, imfilename);
            break;
        }

        // event hadling
        glfwWaitEvents();
    }

    yui::clear_nuklear(nuklear_ctx, legacy_gl);
    yui::clear_glfw(window);
}

int main(int argc, char* argv[]) {
    // command line params
    auto parser = ycmd::make_parser(argc, argv, "view meshes");
    hdr_exposure =
        ycmd::parse_opt<float>(parser, "--exposure", "-e", "image exposure", 0);
    hdr_gamma =
        ycmd::parse_opt<float>(parser, "--gamma", "-g", "image gamma", 2.2);
    amb = ycmd::parse_opt<float>(parser, "--ambient", "", "ambient factor", 0);
    camera_lights = ycmd::parse_flag(parser, "--camera_lights", "-c",
                                     "enable camera lights", false);
    camera = ycmd::parse_opt<int>(parser, "--camera", "-C", "camera", 0);
    no_ui = ycmd::parse_flag(parser, "--no-ui", "", "runs offline", false);
    legacy_gl = ycmd::parse_flag(parser, "--legacy_opengl", "-L",
                                 "uses legacy OpenGL", false);
    aspect = ycmd::parse_opt<float>(parser, "--aspect", "-a", "image aspect",
                                    16.0f / 9.0f);
    res = ycmd::parse_opt<int>(parser, "--resolution", "-r", "image resolution",
                               720);
    imfilename = ycmd::parse_opt<std::string>(parser, "--output", "-o",
                                              "image filename", "out.png");
    filename = ycmd::parse_arg<std::string>(parser, "scene", "scene filename",
                                            "", true);
    ycmd::check_parser(parser);

    // load scene
    scene = yapp::load_scene(filename);
    scene.cameras[camera].aspect = aspect;

    // run ui
    run_ui();

    // done
    return EXIT_SUCCESS;
}
