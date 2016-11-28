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

// view data
struct view_params {
    // filenames
    std::string filename;
    std::string imfilename;

    // window data
    int w = 0, h = 0;
    float exposure = 0, gamma = 2.2f;
    int background = 0;
    ym::vec4f backgrounds[4] = {{0, 0, 0, 0},
                                {0.18f, 0.18f, 0.18f, 0},
                                {0.5f, 0.5f, 0.5f, 0},
                                {1, 1, 1, 0}};

    // drawing
    bool view_edges = false, view_hull = false;

    // mouse state
    int mouse_button = 0;
    ym::vec2f mouse_pos = ym::zero2f;

    // scene
    yapp::scene* scene = nullptr;
    ym::vec3f amb = ym::zero3f;
    int cur_camera = 0;

    // animating
    float time = 0, dt = 1 / 30.0f;
    int frame = 0;
    bool simulating = false;
    ysym::scene* rigid_scene = nullptr;

    // shading
    int shade_prog = 0;
    std::vector<int> shade_txt;

    // lighting
    bool camera_lights;
};

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

ysym::scene make_rigid_scene(yapp::scene& scene, ybvh::scene& scene_bvh) {
    auto rigid_scene = ysym::scene();

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

    return rigid_scene;
}

void simulate_step(yapp::scene& scene, ysym::scene& rigid_scene, float dt) {
    ysym::advance_simulation(rigid_scene, dt);
    for (auto sid = 0; sid < scene.shapes.size(); sid++) {
        scene.shapes[sid].frame = rigid_scene.shapes[sid].frame;
    }
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
    auto aspect = ycmd::parse_opt<float>(parser, "--aspect", "-a",
                                         "image aspect", 16.0f / 9.0f);
    auto dt = ycmd::parse_opt<float>(parser, "--delta_time", "-dt",
                                     "delta time", 1 / 60.0f);
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

    // init rigid simulation
    auto scene_bvh = ybvh::scene();
    auto rigid_scene = make_rigid_scene(scene, scene_bvh);

    // simulation variables
    auto simulating = false;

    // view variables
    auto cur_background = 0;
    const std::array<ym::vec4f, 4> backgrounds = {{{0.0f, 0.0f, 0.0f, 0.0f},
                                                   {0.18f, 0.18f, 0.18f, 0.0f},
                                                   {0.5f, 0.5f, 0.5f, 0.0f},
                                                   {1.0f, 1.0f, 1.0f, 0.0f}}};
    auto wireframe = false;
    auto edges = false;
    auto view_hull = false, view_edges = false;

    // shading
    auto shade_prog = 0;
    std::vector<int> shade_txt;

    // prepare ui context
    auto context = yui::context();

    // init callback
    context.init.push_back(std::function<void(const yui::info& info)>(
        [&](const yui::info& info) {  // load textures
            yapp::init_shade(scene, shade_prog, shade_txt);
        }));

    // window size callback
    context.window_size.push_back(
        std::function<void(const yui::info& info)>([&](const yui::info& info) {
            auto& cam = scene.cameras[camera];
            cam.aspect = (float)info.win_size[0] / (float)info.win_size[1];
        }));

    // window refresh callback
    context.window_refresh.push_back(
        std::function<void(const yui::info& info)>([&](const yui::info& info) {
            // draw
            yapp::shade(scene, camera, shade_prog, shade_txt,
                        backgrounds[cur_background], exposure, gamma, wireframe,
                        edges, camera_lights);
            // draw hull
            if (view_hull) draw_hull(rigid_scene, scene, 1, camera, shade_prog);
        }));

    // check continue callback
    context.update.push_back(
        std::function<int(const yui::info& info)>([&](const yui::info& info) {
            // advance if simulating
            if (simulating) {
                simulate_step(scene, rigid_scene, dt);
                return 1;
            }

            // continue as usual
            return 0;
        }));

    // text callback
    context.text.push_back(
        std::function<void(const yui::info& info, unsigned int)>([&](
            const yui::info& info, unsigned int key) {
            switch (key) {
                case ' ': simulating = !simulating; break;
                case '/': {
                    for (int sid = 0; sid < scene.shapes.size(); sid++) {
                        scene.shapes[sid].frame = ym::identity_frame3f;
                        rigid_scene.shapes[sid].frame = ym::identity_frame3f;
                    }
                } break;
                case '.': simulate_step(scene, rigid_scene, dt); break;
                case '[': exposure -= 1; break;
                case ']': exposure += 1; break;
                case '{': gamma -= 0.1f; break;
                case '}': gamma += 0.1f; break;
                case 'e': view_edges = !view_edges; break;
                case 'h': view_hull = !view_hull; break;
                case 'b':
                    cur_background = (cur_background + 1) % backgrounds.size();
                    break;
                case 's': save_screenshot(info, imfilename); break;
                case 'c': camera_lights = !camera_lights; break;
                case 'C': camera = (camera + 1) % scene.cameras.size(); break;
                default: printf("unsupported key\n"); break;
            }
        }));

    // mouse position callback
    context.mouse_pos.push_back(
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
        }));

    // run ui
    yui::ui_loop(context, (int)std::round(aspect * res), res, "yview");

    // done
    return 0;
}
