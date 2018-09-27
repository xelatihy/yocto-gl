//
// LICENSE:
//
// Copyright (c) 2016 -- 2018 Fabio Pellacini
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

#include "../yocto/ygl.h"
#include "../yocto/yglio.h"
#include "CLI11.hpp"
#include "yglui.h"
using namespace std::literals;

#include <map>

// Application state
struct app_state {
    // scene
    std::shared_ptr<ygl::scene> scn = nullptr;

    // rendering params
    std::string filename = "scene.json"s;
    std::string imfilename = "out.obj"s;
    int camid = 0;                             // camera index
    int resolution = 512;                      // image resolution
    int nsamples = 256;                        // number of samples
    std::string tracer = "pathtrace"s;         // tracer name
    ygl::trace_func tracef = ygl::trace_path;  // tracer
    int nbounces = 4;                          // max depth
    int seed = ygl::trace_default_seed;        // seed
    float pixel_clamp = 100.0f;                // pixel clamping
    int pratio = 8;                            // preview ratio

    // rendering state
    ygl::trace_async_state trace_state = {};

    // view image
    ygl::frame2f imframe = ygl::identity_frame2f;
    bool zoom_to_fit = true;
    float exposure = 0;
    float gamma = 2.2f;
    bool filmic = false;
    ygl::vec4f background = {0.8f, 0.8f, 0.8f, 0};
    ygl::uint gl_txt = 0;
    ygl::uint gl_prog = 0, gl_vbo = 0, gl_ebo;
    bool widgets_open = false;
    ygl::scene_selection selection = {};
    std::vector<ygl::scene_selection> update_list;
    bool navigation_fps = false;
    bool quiet = false;
    int64_t trace_start = 0;
};

auto trace_names = std::vector<std::string>{"pathtrace", "direct",
    "environment", "eyelight", "pathtrace_nomis", "pathtrace_naive",
    "direct_nomis", "debug_normal", "debug_albedo", "debug_texcoord",
    "debug_frontfacing", "debug_diffuse", "debug_specular", "debug_roughness"};

auto tracer_names = std::unordered_map<std::string, ygl::trace_func>{
    {"pathtrace", ygl::trace_path}, {"direct", ygl::trace_direct},
    {"environment", ygl::trace_environment}, {"eyelight", ygl::trace_eyelight},
    {"pathtrace-nomis", ygl::trace_path_nomis},
    {"pathtrace-naive", ygl::trace_path_naive},
    {"direct-nomis", ygl::trace_direct_nomis},
    {"debug_normal", ygl::trace_debug_normal},
    {"debug_albedo", ygl::trace_debug_albedo},
    {"debug_texcoord", ygl::trace_debug_texcoord},
    {"debug_frontfacing", ygl::trace_debug_frontfacing},
    {"debug_diffuse", ygl::trace_debug_diffuse},
    {"debug_specular", ygl::trace_debug_specular},
    {"debug_roughness", ygl::trace_debug_roughness}};

void draw_widgets(GLFWwindow* win, app_state* app) {
    if (ygl::begin_widgets_frame(win, "yitrace", &app->widgets_open)) {
        ImGui::LabelText("scene", "%s", app->filename.c_str());
        ImGui::LabelText("image", "%d x %d @ %d", app->trace_state.img.size.x,
            app->trace_state.img.size.y, app->trace_state.sample);
        if (ImGui::TreeNode("render settings")) {
            auto edited = 0;
            edited += ImGui::Combo("camera", &app->camid,
                app->scn->cameras.size(),
                [&app](int i) { return app->scn->cameras[i]->name.c_str(); });
            edited +=
                ImGui::SliderInt("resolution", &app->resolution, 256, 4096);
            edited += ImGui::SliderInt("nsamples", &app->nsamples, 16, 4096);
            edited += ImGui::Combo("tracer", &app->tracer, trace_names);
            app->tracef = tracer_names.at(app->tracer);
            edited += ImGui::SliderInt("nbounces", &app->nbounces, 1, 10);
            edited += ImGui::SliderInt("seed", (int*)&app->seed, 0, 1000);
            edited += ImGui::SliderInt("pratio", &app->pratio, 1, 64);
            if (edited) app->update_list.push_back(ygl::scene_selection());
            ImGui::LabelText("time/sample", "%0.3lf",
                (app->trace_state.sample) ?
                    (ygl::get_time() - app->trace_start) /
                        (1000000000.0 * app->trace_state.sample) :
                    0.0);
            ImGui::TreePop();
        }
        if (ImGui::TreeNode("view settings")) {
            ImGui::SliderFloat("exposure", &app->exposure, -5, 5);
            ImGui::SliderFloat("gamma", &app->gamma, 1, 3);
            ImGui::ColorEdit4("background", &app->background.x);
            auto zoom = app->imframe.x.x;
            if (ImGui::SliderFloat("zoom", &zoom, 0.1, 10))
                app->imframe.x.x = app->imframe.y.y = zoom;
            ImGui::Checkbox("zoom to fit", &app->zoom_to_fit);
            ImGui::SameLine();
            ImGui::Checkbox("fps", &app->navigation_fps);
            auto mouse_x = 0.0, mouse_y = 0.0;
            glfwGetCursorPos(win, &mouse_x, &mouse_y);
            auto ij = ygl::get_image_coords(
                ygl::vec2f{(float)mouse_x, (float)mouse_y}, app->imframe,
                {app->trace_state.img.size.x, app->trace_state.img.size.y});
            ImGui::DragInt2("mouse", &ij.x);
            if (ij.x >= 0 && ij.x < app->trace_state.img.size.x && ij.y >= 0 &&
                ij.y < app->trace_state.img.size.y) {
                ImGui::ColorEdit4("pixel", &app->trace_state.img[ij].x);
            } else {
                auto zero4f_ = ygl::zero4f;
                ImGui::ColorEdit4("pixel", &zero4f_.x);
            }
            ImGui::TreePop();
        }
        if (ImGui::TreeNode("scene tree")) {
            if (ImGui::Button("print stats")) ygl::print_stats(app->scn);
            ygl::draw_glwidgets_scene_tree(
                "", app->scn, app->selection, app->update_list, 200, {});
            ImGui::TreePop();
        }
        if (ImGui::TreeNode("scene object")) {
            ygl::draw_glwidgets_scene_inspector(
                "", app->scn, app->selection, app->update_list, 200, {});
            ImGui::TreePop();
        }
    }
    ygl::end_widgets_frame();
}

void draw(GLFWwindow* win) {
    auto app = (app_state*)glfwGetWindowUserPointer(win);
    ygl::draw_glimage(win, app->trace_state.display, app->imframe,
        app->zoom_to_fit, app->background);
    draw_widgets(win, app);
    glfwSwapBuffers(win);
}

bool update(const std::shared_ptr<app_state>& app) {
    // exit if no updated
    if (app->update_list.empty()) return false;

    // stop renderer
    ygl::trace_async_stop(app->trace_state);

    // update BVH
    for (auto sel : app->update_list) {
        if (sel.as<ygl::shape>()) {
            ygl::refit_bvh(sel.as<ygl::shape>());
            ygl::refit_bvh(app->scn, false);
        }
        if (sel.as<ygl::instance>()) { ygl::refit_bvh(app->scn, false); }
        if (sel.as<ygl::node>()) {
            ygl::update_transforms(app->scn, 0);
            ygl::refit_bvh(app->scn, false);
        }
    }
    app->update_list.clear();

    app->tracef = tracer_names.at(app->tracer);
    app->trace_state = {};
    app->trace_start = ygl::get_time();
    ygl::trace_async_start(app->trace_state, app->scn, app->camid,
        app->resolution, app->nsamples, app->tracef, app->exposure, app->gamma,
        app->filmic, app->pratio, app->nbounces, app->pixel_clamp, app->seed);

    // updated
    return true;
}

// run ui loop
void run_ui(const std::shared_ptr<app_state>& app) {
    // window
    auto win_size = ygl::clamp(app->trace_state.img.size, 512, 1024);
    auto win = ygl::make_window(win_size, "yitrace", app.get(), draw);

    // init widget
    ygl::init_widgets(win);

    // loop
    auto mouse_pos = ygl::zero2f, last_pos = ygl::zero2f;
    auto mouse_button = 0;
    while (!glfwWindowShouldClose(win)) {
        last_pos = mouse_pos;
        glfwGetCursorPosExt(win, &mouse_pos.x, &mouse_pos.y);
        mouse_button = glfwGetMouseButtonIndexExt(win);
        auto alt_down = glfwGetAltKeyExt(win);
        auto shift_down = glfwGetShiftKeyExt(win);
        auto widgets_active = ImGui::GetWidgetsActiveExt();

        // handle mouse and keyboard for navigation
        if (mouse_button && !alt_down && !widgets_active) {
            auto dolly = 0.0f;
            auto pan = ygl::zero2f;
            auto rotate = ygl::zero2f;
            if (mouse_button == 1) rotate = (mouse_pos - last_pos) / 100.0f;
            if (mouse_button == 2) dolly = (mouse_pos.x - last_pos.x) / 100.0f;
            if (mouse_button == 3 || (mouse_button == 1 && shift_down))
                pan = (mouse_pos - last_pos) / 100.0f;
            auto cam = app->scn->cameras.at(app->camid);
            ygl::camera_turntable(cam->frame, cam->focus, rotate, dolly, pan);
            app->update_list.push_back(cam);
        }

        // selection
        if (mouse_button && alt_down && !widgets_active) {
            auto ij = ygl::get_image_coords(
                mouse_pos, app->imframe, app->trace_state.img.size);
            if (ij.x < 0 || ij.x >= app->trace_state.img.size.x || ij.y < 0 ||
                ij.y >= app->trace_state.img.size.y) {
                auto cam = app->scn->cameras.at(app->camid);
                auto ray = eval_camera_ray(cam, ij,
                    ygl::vec2i{app->trace_state.img.size.x,
                        app->trace_state.img.size.y},
                    {0.5f, 0.5f}, ygl::zero2f);
                auto isec = intersect_ray(app->scn, ray);
                if (isec.ist) app->selection = isec.ist;
            }
        }

        // update
        update(app);

        // draw
        draw(win);

        // event hadling
        if (!mouse_button && !widgets_active)
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        glfwPollEvents();
    }
}

int main(int argc, char* argv[]) {
    // command line parameters
    auto filename = "scene.json"s;  // input filename
    auto imfilename = "out.obj"s;   // image filewname
    auto camid = 0;                 // camera index
    auto resolution = 512;          // image resolution
    auto nsamples = 256;            // number of samples
    auto tracer = "pathtrace"s;     // tracer name
    auto nbounces = 8;              // max depth
    auto seed = 7;                  // seed
    auto pixel_clamp = 100.0f;      // pixel clamping
    auto double_sided = false;      // double sided
    auto add_skyenv = false;        // add sky environment
    auto pratio = 8;                // preview ratio
    auto quiet = false;             // quiet mode

    // parse command line
    CLI::App parser("progressive path tracing", "yitrace");
    parser.add_option("--camera", camid, "Camera index.");
    parser.add_option(
        "--resolution,-r", resolution, "Image vertical resolution.");
    parser.add_option("--nsamples,-s", nsamples, "Number of samples.");
    parser.add_option("--tracer,-t", tracer, "Trace type.")
        ->check([](const std::string& s) -> std::string {
            if (tracer_names.find(s) == tracer_names.end())
                throw CLI::ValidationError("unknown tracer name");
            return s;
        });
    parser.add_option("--nbounces", nbounces, "Maximum number of bounces.");
    parser.add_option("--pixel-clamp", pixel_clamp, "Final pixel clamping.");
    parser.add_option("--seed", seed, "Seed for the random number generators.");
    parser.add_flag(
        "--double-sided,-D", double_sided, "Double-sided rendering.");
    parser.add_flag("--add-skyenv,-E", add_skyenv, "add missing env map");
    parser.add_flag("--quiet,-q", quiet, "Print only errors messages");
    parser.add_option(
        "--pration", pratio, "Preview ratio for async rendering.");
    parser.add_option("--output-image,-o", imfilename, "Image filename");
    parser.add_option("scene", filename, "Scene filename")->required(true);
    try {
        parser.parse(argc, argv);
    } catch (const CLI::ParseError& e) { return parser.exit(e); }

    // scene loading
    auto scn = std::shared_ptr<ygl::scene>();
    if (!quiet) std::cout << "loading scene" << filename << "\n";
    auto load_start = ygl::get_time();
    try {
        scn = ygl::load_scene(filename);
    } catch (const std::exception& e) {
        std::cout << "cannot load scene " << filename << "\n";
        std::cout << "error: " << e.what() << "\n";
        exit(1);
    }
    if (!quiet)
        std::cout << "loading in "
                  << ygl::format_duration(ygl::get_time() - load_start) << "\n";

    // tesselate
    if (!quiet) std::cout << "tesselating scene elements\n";
    ygl::update_tesselation(scn);

    // update bbox and transforms
    ygl::update_transforms(scn);
    ygl::update_bbox(scn);

    // add components
    if (!quiet) std::cout << "adding scene elements\n";
    if (add_skyenv && scn->environments.empty()) {
        scn->environments.push_back(ygl::make_sky_environment("sky"));
        scn->textures.push_back(scn->environments.back()->ke_txt);
    }
    if (double_sided)
        for (auto mat : scn->materials) mat->double_sided = true;
    if (scn->cameras.empty())
        scn->cameras.push_back(ygl::make_bbox_camera("<view>", scn->bbox));
    ygl::add_missing_names(scn);
    for (auto err : ygl::validate(scn)) std::cout << "warning: " << err << "\n";

    // build bvh
    if (!quiet) std::cout << "building bvh\n";
    auto bvh_start = ygl::get_time();
    ygl::update_bvh(scn);
    if (!quiet)
        std::cout << "building bvh in "
                  << ygl::format_duration(ygl::get_time() - bvh_start) << "\n";

    // init renderer
    if (!quiet) std::cout << "initializing lights\n";
    ygl::update_lights(scn);

    // fix renderer type if no lights
    if (scn->lights.empty() && scn->environments.empty() &&
        tracer != "eyelight") {
        if (!quiet)
            std::cout << "no lights presents, switching to eyelight shader\n";
        tracer = "eyelight";
    }

    // prepare application
    auto app = std::make_shared<app_state>();
    app->scn = scn;
    app->filename = filename;
    app->imfilename = imfilename;
    app->camid = camid;
    app->resolution = resolution;
    app->nsamples = nsamples;
    app->tracer = tracer;
    app->tracef = tracer_names.at(tracer);
    app->nbounces = nbounces;
    app->seed = seed;
    app->pixel_clamp = pixel_clamp;
    app->pratio = pratio;

    // initialize rendering objects
    if (!quiet) std::cout << "starting async renderer\n";
    app->trace_start = ygl::get_time();
    ygl::trace_async_start(app->trace_state, app->scn, app->camid,
        app->resolution, app->nsamples, app->tracef, app->exposure, app->gamma,
        app->filmic, app->pratio, app->nbounces, app->pixel_clamp, app->seed);

    // run interactive
    run_ui(app);

    // cleanup
    ygl::trace_async_stop(app->trace_state);

    // done
    return 0;
}
