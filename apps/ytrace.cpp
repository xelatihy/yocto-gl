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
#include "../yocto/yocto_math.h"
#include "../yocto/yocto_shape.h"
#include "../yocto/yocto_trace.h"

#include "ThreadPool.h"

#define MAX_BLOCKS_PER_STEP 4096

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
                const ym::image_view<ym::vec4f>& img) {
    auto ext = ycmd::get_extension(filename);
    if (ext == ".hdr") {
        stbi_write_hdr(filename.c_str(), img.size()[0], img.size()[1], 4,
                       (float*)img.data());
    } else if (ext == ".png") {
        auto img_ub = ym::image<ym::vec4b>(img.size());
        for (auto j = 0; j < img.size()[1]; j++) {
            for (int i = 0; i < img.size()[0]; i++) {
                auto& p = img[{i, j}];
                auto& pub = img_ub[{i, j}];
                pub[0] =
                    ym::clamp((int)(256 * std::pow(p[0], 1 / 2.2f)), 0, 255);
                pub[1] =
                    ym::clamp((int)(256 * std::pow(p[1], 1 / 2.2f)), 0, 255);
                pub[2] =
                    ym::clamp((int)(256 * std::pow(p[2], 1 / 2.2f)), 0, 255);
                pub[3] = ym::clamp((int)(256 * p[3]), 0, 255);
            }
        }
        stbi_write_png(filename.c_str(), img.size()[0], img.size()[1], 4,
                       img_ub.data(), img.size()[0] * 4);
    } else {
        printf("supports only hdr and png for image writing\n");
        return;
    }
}

ytrace::scene init_trace_cb(yapp::scene& scene, ybvh::scene& scene_bvh,
                            int camera) {
    auto trace_scene = ytrace::scene();

    trace_scene.intersect_first = [&scene_bvh](const ym::ray3f& ray) {
        auto isec = ybvh::intersect_first(scene_bvh, ray);
        return reinterpret_cast<ytrace::intersect_point&>(isec);
    };
    trace_scene.intersect_any = [&scene_bvh](const ym::ray3f& ray) {
        return ybvh::intersect_any(scene_bvh, ray);
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

int main(int argc, char* argv[]) {
    auto rtype_names = std::unordered_map<std::string, ytrace::rtype>{
        {"default", ytrace::rtype::def},
        {"uniform", ytrace::rtype::uniform},
        {"stratified", ytrace::rtype::stratified},
        {"cmjs", ytrace::rtype::cmjs}};
    auto stype_names = std::unordered_map<std::string, ytrace::stype>{
        {"default", ytrace::stype::def},
        {"eye", ytrace::stype::eyelight},
        {"direct", ytrace::stype::direct},
        {"path", ytrace::stype::pathtrace}};

    // params
    auto parser = ycmd::make_parser(argc, argv, "trace meshes");
    auto exposure =
        ycmd::parse_opt<float>(parser, "--exposure", "-e", "image exposure", 0);
    auto gamma =
        ycmd::parse_opt<float>(parser, "--gamma", "-g", "image gamma", 2.2);
    auto rtype = ycmd::parse_opte<ytrace::rtype>(
        parser, "--random", "", "random type", ytrace::rtype::def, rtype_names);
    auto stype = ycmd::parse_opte<ytrace::stype>(
        parser, "--integrator", "-i", "integrator type", ytrace::stype::def,
        stype_names);
    auto amb =
        ycmd::parse_opt<float>(parser, "--ambient", "", "ambient factor", 0);
    auto camera_lights = ycmd::parse_flag(parser, "--camera_lights", "-c",
                                          "enable camera lights", false);
    auto camera = ycmd::parse_opt<int>(parser, "--camera", "-C", "camera", 0);
    auto no_ui = ycmd::parse_flag(parser, "--no-ui", "", "runs offline", false);
    auto nthreads = ycmd::parse_opt<int>(
        parser, "--threads", "-t", "number of threads [0 for default]", 0);
    auto block_size =
        ycmd::parse_opt<int>(parser, "--block_size", "", "block size", 32);
    auto samples =
        ycmd::parse_opt<int>(parser, "--samples", "-s", "image samples", 256);
    auto aspect = ycmd::parse_opt<float>(parser, "--aspect", "-a",
                                         "image aspect", 16.0f / 9.0f);
    auto res = ycmd::parse_opt<int>(parser, "--resolution", "-r",
                                    "image resolution", 720);
    auto imfilename = ycmd::parse_opt<std::string>(parser, "--output", "-o",
                                                   "image filename", "out.hdr");
    auto filename = ycmd::parse_arg<std::string>(parser, "scene",
                                                 "scene filename", "", true);
    ycmd::check_parser(parser);

    // setting up multithreading
    if (!nthreads) nthreads = std::thread::hardware_concurrency();
    ThreadPool pool(nthreads);
    std::vector<std::future<void>> futures;

    // loading scene
    auto scene = yapp::load_scene(filename);
    scene.cameras[camera].aspect = aspect;

    // preparing raytracer
    auto scene_bvh = yapp::make_bvh(scene);
    auto trace_scene = init_trace_cb(scene, scene_bvh, camera);
    ytrace::render_params params;
    params.stype = (camera_lights) ? ytrace::stype::eyelight : stype;
    params.rtype = rtype;
    params.amb = {amb, amb, amb};

    // image rendering params
    auto width = (int)std::round(aspect * res), height = res;
    auto blocks = make_image_blocks(width, height, block_size);
    auto img = ym::image<ym::vec4f>({width, height}, ym::zero4f);
    auto buf = ym::image<ym::vec4f>({width, height}, ym::zero4f);

    // launching renderer
    if (no_ui) {
        printf("tracing %s to %s\n", filename.c_str(), imfilename.c_str());
        printf("rendering ...");
        fflush(stdout);
        for (auto cur_sample = 0; cur_sample < samples; cur_sample++) {
            printf("\rrendering sample %d/%d", cur_sample + 1, samples);
            fflush(stdout);
            futures.clear();
            for (auto cur_block = 0; cur_block < blocks.size(); cur_block++) {
                futures.push_back(
                    pool.enqueue([=, &trace_scene, &img, &blocks]() {
                        ytrace::trace_block(
                            trace_scene, camera, img, samples,
                            blocks[cur_block].first, blocks[cur_block].second,
                            {cur_sample, cur_sample + 1}, params, true);
                    }));
            }
            for (auto& future : futures) future.wait();
        }
        printf("\rrendering done\n");
        fflush(stdout);
        // TODO: tonemap
        save_image(imfilename, img);
    } else {
        // prepare ui context
        auto context = yui::context();

        // view data
        auto scene_updated = true;
        auto cur_background = 0;
        const std::array<ym::vec4f, 4> backgrounds = {
            {{0.0f, 0.0f, 0.0f, 0.0f},
             {0.18f, 0.18f, 0.18f, 0.0f},
             {0.5f, 0.5f, 0.5f, 0.0f},
             {1.0f, 1.0f, 1.0f, 0.0f}}};
        auto texture_id = 0;

        // progressize rendering data
        auto preview = ym::image<ym::vec4f>(
            {width / block_size, height / block_size}, ym::zero4f);
        auto cur_sample = 0, cur_block = 0;
        auto blocks_per_update = 8;

        // init callback
        context.init.push_back(std::function<void(const yui::info& info)>(
            [&](const yui::info& info) {
                texture_id = yglu::make_texture(
                    width, height, 4, (float*)img.data(), false, true);
            }));

        // window size callback
        context.window_size.push_back(
            std::function<void(const yui::info& info)>(
                [&](const yui::info& info) {
                    auto& cam = scene.cameras[camera];
                    width = info.win_size[0];
                    height = info.win_size[1];
                    cam.aspect = (float)width / (float)height;
                    // updated images
                    img.resize({width, height});
                    buf.resize({width, height});
                    preview.resize({width / block_size, height / block_size});

                    // update texture
                    if (texture_id) yglu::clear_texture(&texture_id);
                    texture_id = yglu::make_texture(
                        width, height, 4, (float*)img.data(), false, true);

                    // updated blocks
                    blocks = make_image_blocks(width, height, block_size);
                    scene_updated = true;
                }));

        // window refresh callback
        context.window_refresh.push_back(
            std::function<void(const yui::info& info)>(
                [&](const yui::info& info) {
                    //                char title[4096];
                    //                sprintf(title, "ytrace | %dx%d | %d/%d  |
                    //                %s",
                    //                width, height, cur_sample,
                    //                        view->ns, view->filename.c_str());
                    //                glfwSetWindowTitle(window, title);

                    auto background = backgrounds[cur_background];
                    glClearColor(background[0], background[1], background[2],
                                 background[3]);
                    glDisable(GL_DEPTH_TEST);
                    glClear(GL_COLOR_BUFFER_BIT);

                    yglu::shade_image(texture_id, width, height, width, height,
                                      0, 0, 1, exposure, gamma);
                }));

        // check continue callback
        context.update.push_back(std::function<int(const yui::info& info)>([&](
            const yui::info& info) {
            if (scene_updated) {
                // update camera
                auto& cam = scene.cameras[camera];
                auto& trace_cam = trace_scene.cameras[camera];
                trace_cam = {cam.frame, cam.yfov, cam.aspect, cam.aperture,
                             cam.focus};

                // render preview
                ytrace::trace_image(trace_scene, camera, preview, 1, params);
                for (auto qj = 0; qj < height / block_size; qj++) {
                    for (auto qi = 0; qi < width / block_size; qi++) {
                        for (auto j = qj * block_size;
                             j < ym::min((qj + 1) * block_size, height); j++) {
                            for (auto i = qi * block_size;
                                 i < ym::min((qi + 1) * block_size, width);
                                 i++) {
                                img[{i, j}] = preview[{qi, qj}];
                            }
                        }
                    }
                }
                yglu::update_texture(texture_id, width, height, 4,
                                     (float*)img.data());

                // reset current counters
                cur_sample = 0;
                cur_block = 0;
                scene_updated = false;
            } else {
                if (cur_sample == samples) return 0;
                futures.clear();
                for (auto b = 0;
                     cur_block < blocks.size() && b < blocks_per_update;
                     cur_block++, b++) {
                    futures.push_back(pool.enqueue([=, &trace_scene, &img,
                                                    &blocks]() {
                        ytrace::trace_block(
                            trace_scene, camera, img, samples,
                            blocks[cur_block].first, blocks[cur_block].second,
                            {cur_sample, cur_sample + 1}, params, true);
                    }));
                }
                for (auto& future : futures) future.wait();
                yglu::update_texture(texture_id, width, height, 4,
                                     (float*)img.data());
                if (cur_block == blocks.size()) {
                    cur_block = 0;
                    cur_sample++;
                }
            }

            return 1;
        }));

        // text callback
        context.text.push_back(
            std::function<void(const yui::info& info, unsigned int)>(
                [&](const yui::info& info, unsigned int key) {
                    switch (key) {
                        case '1':
                            exposure = 1;
                            gamma = 1;
                            break;
                        case '2':
                            exposure = 1;
                            gamma = 2.2f;
                            break;
                        case '[': exposure -= 1; break;
                        case ']': exposure += 1; break;
                        case '{': gamma -= 0.1f; break;
                        case '}': gamma += 0.1f; break;
                        case 's': {
                            save_image(imfilename.c_str(), img);
                        } break;
                        case 'C': {
                            camera = (camera + 1) % scene.cameras.size();
                            scene_updated = true;
                        } break;
                        case 'b':
                            cur_background =
                                (cur_background + 1) % backgrounds.size();
                            break;
                        default: printf("unsupported key\n"); break;
                    }
                }));

        // mouse position callback
        context.mouse_pos.push_back(std::function<void(const yui::info& info)>(
            [&](const yui::info& info) {
                if (info.mouse_button) {
                    auto dolly = 0.0f;
                    auto pan = ym::zero2f;
                    auto rotate = ym::zero2f;
                    switch (info.mouse_button) {
                        case 1:
                            rotate = (info.mouse_pos - info.mouse_last) / 100;
                            break;
                        case 2:
                            dolly = (info.mouse_pos[0] - info.mouse_last[0]) /
                                    100.0f;
                            break;
                        case 3:
                            pan = (info.mouse_pos - info.mouse_last) / 100;
                            break;
                        default: break;
                    }

                    auto& cam = scene.cameras[camera];
                    ym::turntable(cam.frame, cam.focus, rotate, dolly, pan);

                    scene_updated = true;
                }
            }));

        // run ui
        yui::ui_loop(context, (int)std::round(aspect * res), res, "yview");
    }

    // done
    return EXIT_SUCCESS;
}
