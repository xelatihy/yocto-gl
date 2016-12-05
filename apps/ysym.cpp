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
#include "../yocto/yocto_sym.h"

// scene
std::string filename;
std::string imfilename;
yapp::scene scene;

// simulating
ysym::scene rigid_scene;
ybvh::scene scene_bvh;
bool simulating = false;
float dt = 1 / 60.f;
int frame = 0;

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
    // draw hull
    // if (hull) draw_hull(rigid_scene, scene, 1, camera, shade_prog);
}

// deprecated function removed from build
#if 0
void draw_hull(ysym::scene* rigid_scene, yapp::scene* scene, float dt,
               int cur_camera, int prog) {
    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    yo_camera* cam = &scene->cameras[cur_camera];
    yocto::ym::frame3f camera_frame = yocto::ym::lookat_xform3(
        *(yocto::ym::vec3f*)cam->from, *(yocto::ym::vec3f*)cam->to, *(yocto::ym::vec3f*)cam->up);
    yocto::ym::mat4f camera_xform = yocto::ym::mat4f(camera_frame);
    yocto::ym::mat4f camera_view = yocto::ym::mat4f(yocto::ym::inverse(camera_frame));
    yocto::ym::mat4f camera_proj = yocto::ym::perspective_mat4(
        2 * atanf(cam->height / 2), cam->width / cam->height, 0.1f, 10000.0f);

    yglu::stdshader_begin_frame(prog, true, 1, 2.2, camera_xform.data(),
                             camera_view.data(), camera_proj.data());

    for (int i = 0; i < scene->nshapes; i++) {
        yo_shape* shape = &scene->shapes[i];
        ysym::_body* body = &rigid_scene->bodies[i];

        yglu::stdshader_begin_shape(prog, shape->xform);

        float zero[3] = {0, 0, 0};
        float one[3] = {1, 1, 1};
        yglu::stdshader_set_material(prog, shape->etype, one, zero, zero, 0, 0, 0,
                                  0, 0, false);

        glLineWidth(2);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        yglu::stdshader_set_vert(prog, shape->pos, 0, 0, 0);
        yglu::stdshader_draw_elem(prog, body->shape.nelems, (int*)body->shape.elem,
                               body->shape.etype);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glLineWidth(1);

        yglu::stdshader_end_shape();
    }

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

        yocto::ym::vec3f pos[] = {col->frame.pos,
                          col->frame.pos + col->depth * col->frame.norm};
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

    yglu::stdshader_end_frame();

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glEnable(GL_DEPTH_TEST);
}
#else
void draw_hull(const ysym::scene& rigid_scene, const yapp::scene& scene,
               float dt, int cur_camera, int prog) {}
#endif

void make_rigid_scene(yapp::scene& scene, ysym::scene& rigid_scene,
                      ybvh::scene& scene_bvh) {
    // add each shape
    auto first_hack = true;
    for (auto& shape : scene.shapes) {
        auto& mat = scene.materials[shape.matid];
        auto simulated =
            !first_hack && ym::length(mat.ke) == 0 && !shape.triangles.empty();
        auto density = (simulated) ? 1.0f : 0.0f;
        rigid_scene.shapes.push_back({shape.frame, ym::zero3f, ym::zero3f,
                                      density, simulated, shape.triangles,
                                      shape.pos});
        first_hack = false;
    }

    // set up final bvh
    scene_bvh = ybvh::scene();
    for (auto i = 0; i < scene.shapes.size(); i++) {
        auto& shape = scene.shapes[i];
        assert(!shape.points.empty() || !shape.lines.empty() ||
               !shape.triangles.empty());
        scene_bvh.shapes.push_back({ym::to_mat(shape.frame),
                                    ym::to_mat(ym::inverse(shape.frame)),
                                    shape.points, shape.lines, shape.triangles,
                                    shape.pos, shape.radius});
    }
    ybvh::build_bvh(scene_bvh);

    // setup collisions
    rigid_scene.overlap_shapes =
        [&scene_bvh](std::vector<ym::vec2i>& overlaps) {
            return ybvh::overlap_shape_bounds(scene_bvh, true, overlaps);
        };
    rigid_scene.overlap_shape = [&scene_bvh](int sid, const ym::vec3f& pt,
                                             float max_dist) {
        auto overlap = ybvh::overlap_first(scene_bvh.shapes[sid], pt, max_dist);
        return *(ysym::overlap_point*)&overlap;
    };
    rigid_scene.overlap_refit = [&scene_bvh, &rigid_scene]() {
        for (auto sid = 0; sid < rigid_scene.shapes.size(); sid++) {
            scene_bvh.shapes[sid].xform =
                ym::to_mat(rigid_scene.shapes[sid].frame);
            scene_bvh.shapes[sid].inv_xform =
                ym::to_mat(ym::inverse(rigid_scene.shapes[sid].frame));
        }
        ybvh::refit_bvh(scene_bvh);
    };

    // initialize
    ysym::init_simulation(rigid_scene);
}

void simulate_step(yapp::scene& scene, ysym::scene& rigid_scene, float dt) {
    ysym::advance_simulation(rigid_scene, dt);
    for (auto sid = 0; sid < scene.shapes.size(); sid++) {
        scene.shapes[sid].frame = rigid_scene.shapes[sid].frame;
    }
}

void text_callback(GLFWwindow* window, unsigned int key) {
    nk_glfw3_gl3_char_callback(window, key);
    if (nk_item_is_any_active(nuklear_ctx)) return;
    switch (key) {
        case ' ': simulating = !simulating; break;
        case '/': {
            for (int sid = 0; sid < scene.shapes.size(); sid++) {
                scene.shapes[sid].frame = ym::identity_frame3f;
                rigid_scene.shapes[sid].frame = ym::identity_frame3f;
            }
            frame = 0;
        } break;
        case '.':
            simulate_step(scene, rigid_scene, dt);
            frame += 1;
            break;
        case '[': hdr_exposure -= 1; break;
        case ']': hdr_exposure += 1; break;
        case '{': hdr_gamma -= 0.1f; break;
        case '}': hdr_gamma += 0.1f; break;
        case 'e': edges = !edges; break;
        case 'w': wireframe = !wireframe; break;
        case 'b':
            cur_background = (cur_background + 1) % backgrounds.size();
            break;
        case 's': save_screenshot(window, imfilename); break;
        case 'c': camera_lights = !camera_lights; break;
        case 'C': camera = (camera + 1) % scene.cameras.size(); break;
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
    if (nk_begin(nuklear_ctx, "ysym", nk_rect(window_size[0] - hud_width, 0,
                                              hud_width, window_size[1]),
                 NK_WINDOW_BORDER | NK_WINDOW_MOVABLE | NK_WINDOW_SCALABLE |
                     NK_WINDOW_MINIMIZABLE | NK_WINDOW_TITLE)) {
        nk_layout_row_dynamic(nuklear_ctx, 30, 1);
        nk_label(nuklear_ctx, filename.c_str(), NK_TEXT_LEFT);
        nk_layout_row_dynamic(nuklear_ctx, 30, 2);
        nk_value_int(nuklear_ctx, "frame", frame);
        nk_property_float(nuklear_ctx, "dt", 0, &dt, 1, 1 / 240.0f, 1 / 240.0f);
        if (nk_button_label(nuklear_ctx, "start")) simulating = true;
        if (nk_button_label(nuklear_ctx, "stop")) simulating = false;
        if (nk_button_label(nuklear_ctx, "step")) {
            simulate_step(scene, rigid_scene, dt);
            frame += 1;
        }
        if (nk_button_label(nuklear_ctx, "reset")) {
            for (int sid = 0; sid < scene.shapes.size(); sid++) {
                scene.shapes[sid].frame = ym::identity_frame3f;
                rigid_scene.shapes[sid].frame = ym::identity_frame3f;
            }
            frame = 0;
        }
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
    window = yui::init_glfw({(int)(aspect * res), res}, "ysym", legacy_gl,
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

        // advance if simulating
        if (simulating) {
            simulate_step(scene, rigid_scene, dt);
            frame += 1;
        }

        // check for screenshot
        if (no_ui) {
            save_screenshot(window, imfilename);
            break;
        }

        // event hadling
        if (simulating)
            glfwPollEvents();
        else
            glfwWaitEvents();
    }

    yui::clear_nuklear(nuklear_ctx, legacy_gl);
    yui::clear_glfw(window);
}

int main(int argc, char* argv[]) {
    // command line
    auto parser = ycmd::make_parser(argc, argv, "view meshes");
    hdr_exposure =
        ycmd::parse_opt<float>(parser, "--exposure", "-e", "image exposure", 0);
    hdr_gamma =
        ycmd::parse_opt<float>(parser, "--gamma", "-g", "image gamma", 2.2);
    amb = ycmd::parse_opt<float>(parser, "--ambient", "", "ambient factor", 0);
    camera_lights = ycmd::parse_flag(parser, "--camera_lights", "-c",
                                     "enable camera lights", false);
    camera = ycmd::parse_opt<int>(parser, "--camera", "-C", "camera", 0);
    legacy_gl = ycmd::parse_flag(parser, "--legacy_opengl", "-L",
                                 "uses legacy OpenGL", false);
    aspect = ycmd::parse_opt<float>(parser, "--aspect", "-a", "image aspect",
                                    16.0f / 9.0f);
    dt = ycmd::parse_opt<float>(parser, "--delta_time", "-dt", "delta time",
                                1 / 60.0f);
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

    // init rigid simulation
    make_rigid_scene(scene, rigid_scene, scene_bvh);

    // run ui
    run_ui();

    // done
    return 0;
}
