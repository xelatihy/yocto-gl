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

#include "../yocto/yocto_bvh.h"
#include "../yocto/yocto_cmd.h"
#include "../yocto/yocto_trace.h"

#include "ThreadPool.h"

// scene
std::string filename;
std::string imfilename;
yapp::scene scene;

// rendered image
ym::image<ym::vec4f> hdr;
ym::image<ym::vec4b> ldr;

// rendering scene
ybvh::scene scene_bvh;
ytrace::scene trace_scene;

// rendering parameters
int samples = 256;
ytrace::render_params params;
int block_size = 32;
int nthreads = 0;

// lighting
float hdr_exposure = 0;
float hdr_gamma = 2.2;
bool camera_lights = false;

// camera
int camera = 0;
float aspect = 16.0f / 9.0f;
int res = 720;

// gl
bool no_ui = false;
bool legacy_gl = false;

// progressive rendering
ym::image<ym::vec4f> preview;
int cur_sample = 0, cur_block = 0;
int blocks_per_update = 8;

// view variables
int cur_background = 0;
const std::array<ym::vec4f, 4> backgrounds = {{{0.0f, 0.0f, 0.0f, 0.0f},
                                               {0.18f, 0.18f, 0.18f, 0.0f},
                                               {0.5f, 0.5f, 0.5f, 0.0f},
                                               {1.0f, 1.0f, 1.0f, 0.0f}}};
bool wireframe = false;
bool edges = false;

// shading
bool scene_updated = true;
yglu::uint texture_id = 0;
float texture_exposure, texture_gamma;

// cacahed rendering values
std::vector<std::pair<ym::vec2i, ym::vec2i>> blocks;
ThreadPool* pool = nullptr;
std::vector<std::future<void>> futures;

// widget
GLFWwindow* window = nullptr;

// nuklear
nk_context* nuklear_ctx = nullptr;
int hud_width = 256;

std::vector<std::pair<ym::vec2i, ym::vec2i>> make_image_blocks(int w, int h,
                                                               int bs) {
    std::vector<std::pair<ym::vec2i, ym::vec2i>> blocks;
    for (int j = 0; j < h; j += bs) {
        for (int i = 0; i < w; i += bs) {
            blocks.push_back(
                {{i, j}, {ym::min(bs, w - i), ym::min(bs, h - j)}});
        }
    }
    return blocks;
}

void save_image(const std::string& filename,
                const ym::image_view<ym::vec4f>& hdr,
                const ym::image_view<ym::vec4b>& ldr) {
    auto ext = ycmd::get_extension(filename);
    if (ext == ".hdr") {
        stbi_write_hdr(filename.c_str(), hdr.size()[0], hdr.size()[1], 4,
                       (float*)hdr.data());
    } else if (ext == ".png") {
        stbi_write_png(filename.c_str(), ldr.size()[0], ldr.size()[1], 4,
                       ldr.data(), ldr.size()[0] * 4);
    } else {
        printf("supports only hdr and png for image writing\n");
        return;
    }
}

ybvh::scene make_bvh(const yapp::scene& scene) {
    auto scene_bvh = ybvh::scene();
    auto sid = 0;
    for (auto& shape : scene.shapes) {
        scene_bvh.shapes.push_back({sid++,
                                    shape.frame,
                                    shape.points,
                                    shape.lines,
                                    shape.triangles,
                                    {},
                                    shape.pos,
                                    shape.radius});
    }
    ybvh::build_bvh(scene_bvh);
    return scene_bvh;
}

ytrace::scene make_trace_scene(yapp::scene& scene, ybvh::scene& scene_bvh,
                               int camera) {
    auto trace_scene = ytrace::scene();

    trace_scene.intersect_first = [&scene_bvh](const ym::ray3f& ray) {
        auto isec = ybvh::intersect_ray(scene_bvh, ray, false);
        return ytrace::intersect_point{isec.dist,
                                       isec.sid,
                                       isec.eid,
                                       {isec.euv[0], isec.euv[1], isec.euv[2]}};
    };
    trace_scene.intersect_any = [&scene_bvh](const ym::ray3f& ray) {
        return ybvh::intersect_ray(scene_bvh, ray, false);
    };

    for (auto& cam : scene.cameras) {
        trace_scene.cameras.push_back(
            {cam.frame, cam.yfov, cam.aspect, cam.aperture, cam.focus});
    }

    for (auto& env : scene.environments) {
        auto& mat = scene.materials[env.matid];
        trace_scene.environments.push_back({env.frame, mat.ke, mat.ke_txt});
    }

    for (auto& shape : scene.shapes) {
        trace_scene.shapes.push_back(
            {shape.frame,
             shape.matid,
             {shape.points.size(), (ym::vec1i*)shape.points.data()},
             shape.lines,
             shape.triangles,
             shape.pos,
             shape.norm,
             shape.texcoord,
             shape.color,
             shape.radius});
    }

    for (auto& mat : scene.materials) {
        trace_scene.materials.push_back({mat.ke, mat.kd, mat.ks, mat.rs,
                                         ym::zero3f, ym::zero3f, mat.ke_txt,
                                         mat.kd_txt, mat.ks_txt, -1, false});
    }

    for (auto& txt : scene.textures) {
        if (!txt.hdr.empty()) {
            trace_scene.textures.push_back(
                {ym::image_view<ym::vec4f>(txt.hdr), {}});
        } else if (!txt.ldr.empty()) {
            trace_scene.textures.push_back(
                {{}, ym::image_view<ym::vec4b>{txt.ldr}});
        } else
            assert(false);
    }

    ytrace::init_lights(trace_scene);

    return trace_scene;
}

void text_callback(GLFWwindow* window, unsigned int key) {
    nk_glfw3_gl3_char_callback(window, key);
    if (nk_item_is_any_active(nuklear_ctx)) return;
    switch (key) {
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
        case 'b':
            cur_background = (cur_background + 1) % backgrounds.size();
            break;
        case 's': save_image(imfilename.c_str(), hdr, ldr); break;
        default: printf("unsupported key\n"); break;
    }
}

void draw_image() {
    auto framebuffer_size = yui::framebuffer_size(window);
    glViewport(0, 0, framebuffer_size[0], framebuffer_size[1]);

    // begin frame
    auto background = backgrounds[cur_background];
    glClearColor(background[0], background[1], background[2], background[3]);
    glDisable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT);

    // draw image
    auto window_size = yui::window_size(window);
    if (legacy_gl) {
        yglu::legacy::draw_image(texture_id, ldr.size()[0], ldr.size()[1],
                                 window_size[0], window_size[1], 0, 0, 1);
    } else {
        yglu::modern::shade_image(texture_id, ldr.size()[0], ldr.size()[1],
                                  window_size[0], window_size[1], 0, 0, 1);
    }
}

void draw_widgets() {
    auto window_size = yui::window_size(window);
    if (legacy_gl) {
        nk_glfw3_gl2_new_frame();
    } else {
        nk_glfw3_gl3_new_frame();
    }
    if (nk_begin(nuklear_ctx, "ytrace", nk_rect(window_size[0] - hud_width, 0,
                                                hud_width, window_size[1]),
                 NK_WINDOW_BORDER | NK_WINDOW_MOVABLE | NK_WINDOW_SCALABLE |
                     NK_WINDOW_MINIMIZABLE | NK_WINDOW_TITLE)) {
        nk_layout_row_dynamic(nuklear_ctx, 30, 1);
        nk_label(nuklear_ctx, filename.c_str(), NK_TEXT_LEFT);
        nk_layout_row_dynamic(nuklear_ctx, 30, 3);
        nk_value_int(nuklear_ctx, "w", hdr.size()[0]);
        nk_value_int(nuklear_ctx, "h", hdr.size()[1]);
        nk_value_int(nuklear_ctx, "s", cur_sample);
        nk_layout_row_dynamic(nuklear_ctx, 30, 1);
        nk_property_int(nuklear_ctx, "samples", 0, &samples, 1000000, 1, 1);
        nk_layout_row_dynamic(nuklear_ctx, 30, 2);
        nk_property_int(nuklear_ctx, "camera", 0, &camera,
                        (int)scene.cameras.size() - 1, 1, 1);
        camera_lights = nk_check_label(nuklear_ctx, "eyelight", camera_lights);
        nk_layout_row_dynamic(nuklear_ctx, 30, 1);
        nk_property_float(nuklear_ctx, "exposure", -20, &hdr_exposure, 20, 1,
                          1);
        nk_property_float(nuklear_ctx, "gamma", 0.1, &hdr_gamma, 5, 0.1, 0.1);
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

bool update() {
    if (scene_updated) {
        // update camera
        auto& cam = scene.cameras[camera];
        auto& trace_cam = trace_scene.cameras[camera];
        trace_cam = {cam.frame, cam.yfov, cam.aspect, cam.aperture, cam.focus};

        // render preview
        preview.resize(hdr.size() / block_size);
        ytrace::trace_image(trace_scene, camera, preview, 1, params);
        for (auto qj = 0; qj < preview.size()[1]; qj++) {
            for (auto qi = 0; qi < preview.size()[0]; qi++) {
                for (auto j = qj * block_size;
                     j < ym::min((qj + 1) * block_size, hdr.size()[1]); j++) {
                    for (auto i = qi * block_size;
                         i < ym::min((qi + 1) * block_size, hdr.size()[0]);
                         i++) {
                        hdr[{i, j}] = preview[{qi, qj}];
                    }
                }
            }
        }
        ym::exposure_gamma(hdr, ldr, hdr_exposure, hdr_gamma);
        if (legacy_gl) {
            yglu::legacy::update_texture(texture_id, hdr.size()[0],
                                         hdr.size()[1], 4,
                                         (unsigned char*)ldr.data(), false);
        } else {
            yglu::modern::update_texture(texture_id, hdr.size()[0],
                                         hdr.size()[1], 4,
                                         (unsigned char*)ldr.data(), false);
        }

        // reset current counters
        cur_sample = 0;
        cur_block = 0;
        scene_updated = false;
    } else {
        if (cur_sample == samples) return false;
        futures.clear();
        for (auto b = 0; cur_block < blocks.size() && b < blocks_per_update;
             cur_block++, b++) {
            auto block = blocks[cur_block];
            futures.push_back(pool->enqueue([block]() {
                ytrace::trace_block(trace_scene, camera, hdr, samples,
                                    block.first, block.second,
                                    {cur_sample, cur_sample + 1}, params, true);
                ym::exposure_gamma(hdr, ldr, hdr_exposure, hdr_gamma,
                                   block.first, block.second);
            }));
        }
        for (auto& future : futures) future.wait();
        if (texture_exposure != hdr_exposure || texture_gamma != hdr_gamma) {
            ym::exposure_gamma(hdr, ldr, hdr_exposure, hdr_gamma);
            texture_exposure = hdr_exposure;
            texture_gamma = hdr_gamma;
        }
        if (legacy_gl) {
            yglu::legacy::update_texture(texture_id, ldr.size()[0],
                                         ldr.size()[1], 4,
                                         (unsigned char*)ldr.data(), false);
        } else {
            yglu::modern::update_texture(texture_id, ldr.size()[0],
                                         ldr.size()[1], 4,
                                         (unsigned char*)ldr.data(), false);
        }
        if (cur_block == blocks.size()) {
            cur_block = 0;
            cur_sample++;
        }
    }

    return true;
}

void run_ui() {
    // window
    window = yui::init_glfw({(int)(res * aspect), res}, "ytrace", legacy_gl,
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

    nuklear_ctx = yui::init_nuklear(window, legacy_gl);

    if (legacy_gl) {
        texture_id = yglu::legacy::make_texture(hdr.size()[0], hdr.size()[1], 4,
                                                (unsigned char*)ldr.data(),
                                                false, false);
    } else {
        texture_id = yglu::modern::make_texture(hdr.size()[0], hdr.size()[1], 4,
                                                (unsigned char*)ldr.data(),
                                                false, false, false);
    }

    while (!glfwWindowShouldClose(window)) {
        glfwGetWindowSize(window, &window_size[0], &window_size[1]);
        glfwGetFramebufferSize(window, &framebuffer_size[0],
                               &framebuffer_size[1]);

        mouse_last = mouse_pos;
        mouse_pos = yui::mouse_pos(window);
        mouse_button = yui::mouse_button(window);

        // check for window size
        auto window_size = yui::window_size(window);
        if (window_size != hdr.size()) {
            auto& cam = scene.cameras[camera];
            cam.aspect = (float)window_size[0] / (float)window_size[1];
            // updated images
            hdr.resize(window_size);
            ldr.resize(window_size);
            preview.resize(
                {window_size[0] / block_size, window_size[1] / block_size});

            // update texture
            if (legacy_gl) {
                if (texture_id) yglu::legacy::clear_texture(&texture_id);
                texture_id = yglu::legacy::make_texture(
                    hdr.size()[0], hdr.size()[1], 4, (unsigned char*)ldr.data(),
                    false, false);
            } else {
                if (texture_id) yglu::legacy::clear_texture(&texture_id);
                texture_id = yglu::modern::make_texture(
                    hdr.size()[0], hdr.size()[1], 4, (unsigned char*)ldr.data(),
                    false, false, false);
            }

            // updated blocks
            blocks =
                make_image_blocks(hdr.size()[0], hdr.size()[1], block_size);
            scene_updated = true;
        }

        glfwSetWindowTitle(window, ("ytrace | " + filename + " | " +
                                    std::to_string(hdr.size()[0]) + "x" +
                                    std::to_string(hdr.size()[1]) + "@" +
                                    std::to_string(cur_sample))
                                       .c_str());

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
            scene_updated = true;
        }

        // update
        auto updated = update();

        // draw
        draw_image();

        // make ui
        draw_widgets();

        // swap buffers
        glfwSwapBuffers(window);

        // event hadling
        if (updated)
            glfwPollEvents();
        else
            glfwWaitEvents();
    }

    yui::clear_nuklear(nuklear_ctx, legacy_gl);
    yui::clear_glfw(window);
}

void render_offline() {
    printf("tracing %s to %s\n", filename.c_str(), imfilename.c_str());
    printf("rendering ...");
    fflush(stdout);
    for (auto cur_sample = 0; cur_sample < samples; cur_sample++) {
        printf("\rrendering sample %d/%d", cur_sample + 1, samples);
        fflush(stdout);
        futures.clear();
        for (auto cur_block = 0; cur_block < blocks.size(); cur_block++) {
            auto block = blocks[cur_block];
            futures.push_back(pool->enqueue([=, &block]() {
                ytrace::trace_block(trace_scene, camera, hdr, samples,
                                    block.first, block.second,
                                    {cur_sample, cur_sample + 1}, params, true);
            }));
        }
        for (auto& future : futures) future.wait();
    }
    printf("\rrendering done\n");
    fflush(stdout);
    ym::exposure_gamma(hdr, ldr, hdr_exposure, hdr_gamma);
    save_image(imfilename, hdr, ldr);
}

int main(int argc, char* argv[]) {
    auto rtype_names = std::unordered_map<std::string, ytrace::rng_type>{
        {"default", ytrace::rng_type::def},
        {"uniform", ytrace::rng_type::uniform},
        {"stratified", ytrace::rng_type::stratified},
        {"cmjs", ytrace::rng_type::cmjs}};
    auto stype_names = std::unordered_map<std::string, ytrace::shader_type>{
        {"default", ytrace::shader_type::def},
        {"eye", ytrace::shader_type::eyelight},
        {"direct", ytrace::shader_type::direct},
        {"path", ytrace::shader_type::pathtrace}};

    // params
    auto parser = ycmd::make_parser(argc, argv, "trace meshes");
    hdr_exposure =
        ycmd::parse_opt<float>(parser, "--exposure", "-e", "image exposure", 0);
    hdr_gamma =
        ycmd::parse_opt<float>(parser, "--gamma", "-g", "image gamma", 2.2);
    params.rtype = ycmd::parse_opte<ytrace::rng_type>(
        parser, "--random", "", "random type", ytrace::rng_type::def,
        rtype_names);
    params.stype = ycmd::parse_opte<ytrace::shader_type>(
        parser, "--integrator", "-i", "integrator type",
        ytrace::shader_type::def, stype_names);
    camera_lights = ycmd::parse_flag(parser, "--camera_lights", "-c",
                                     "enable camera lights", false);
    camera = ycmd::parse_opt<int>(parser, "--camera", "-C", "camera", 0);
    no_ui = ycmd::parse_flag(parser, "--no-ui", "", "runs offline", false);
    legacy_gl = ycmd::parse_flag(parser, "--legacy_opengl", "-L",
                                 "uses legacy OpenGL", false);
    nthreads = ycmd::parse_opt<int>(parser, "--threads", "-t",
                                    "number of threads [0 for default]", 0);
    block_size =
        ycmd::parse_opt<int>(parser, "--block_size", "", "block size", 32);
    samples =
        ycmd::parse_opt<int>(parser, "--samples", "-s", "image samples", 256);
    aspect = ycmd::parse_opt<float>(parser, "--aspect", "-a", "image aspect",
                                    16.0f / 9.0f);
    res = ycmd::parse_opt<int>(parser, "--resolution", "-r", "image resolution",
                               720);
    imfilename = ycmd::parse_opt<std::string>(parser, "--output", "-o",
                                              "image filename", "out.hdr");
    filename = ycmd::parse_arg<std::string>(parser, "scene", "scene filename",
                                            "", true);
    ycmd::check_parser(parser);

    // setting up multithreading
    if (!nthreads) nthreads = std::thread::hardware_concurrency();
    pool = new ThreadPool(nthreads);

    // loading scene
    scene = yapp::load_scene(filename);
    scene.cameras[camera].aspect = aspect;

    // preparing raytracer
    scene_bvh = make_bvh(scene);
    trace_scene = make_trace_scene(scene, scene_bvh, camera);
    params.stype =
        (camera_lights) ? ytrace::shader_type::eyelight : params.stype;

    // image rendering params
    auto width = (int)std::round(aspect * res), height = res;
    hdr = ym::image<ym::vec4f>({width, height}, ym::zero4f);
    ldr = ym::image<ym::vec4b>({width, height});
    blocks = make_image_blocks(width, height, block_size);
    texture_exposure = hdr_exposure;
    texture_gamma = hdr_gamma;

    // launching renderer
    if (no_ui) {
        render_offline();
    } else {
        // run ui
        run_ui();
    }

    // done
    return EXIT_SUCCESS;
}
