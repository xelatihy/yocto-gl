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
#include "yglutils.h"
#include "ysceneui.h"

// Application state
struct app_state {
    // scene
    ygl::scene* scn = nullptr;

    // rendering params
    std::string filename = "scene.json";
    std::string imfilename = "out.obj";
    int camid = 0;                       // camera id
    int resolution = 512;                // image resolution
    int nsamples = 256;                  // number of samples
    int tracer_id = 0;                   // tracer index
    int nbounces = 4;                    // max depth
    int seed = ygl::trace_default_seed;  // seed
    float pixel_clamp = 100.0f;          // pixel clamping
    int pratio = 8;                      // preview ratio

    // rendering state
    ygl::image4f img = {};
    ygl::image4f display = {};
    std::vector<ygl::rng_state> rng = {};
    bool stop = false;
    std::vector<std::thread> threads;
    int sample = 0;

    // view image
    ygl::vec2f imcenter = ygl::zero2f;
    float imscale = 1;
    bool zoom_to_fit = true;
    float exposure = 0;
    float gamma = 2.2f;
    bool filmic = false;
    ygl::vec4f background = {0.8f, 0.8f, 0.8f, 0};
    unsigned int gl_txt = 0;
    unsigned int gl_prog = 0, gl_vbo = 0, gl_ebo;
    bool widgets_open = false;
    void* selection = nullptr;
    std::vector<std::pair<std::string, void*>> update_list;
    bool navigation_fps = false;
    bool quiet = false;
    int64_t trace_start = 0;

    ~app_state() {
        if (scn) delete scn;
    }
};

auto tracer_names = std::vector<std::string>{"pathtrace", "pathtrace_volume",
    "direct", "environment", "eyelight", "pathtrace-nomis", "pathtrace-naive",
    "direct-nomis", "debug_normal", "debug_albedo", "debug_texcoord",
    "debug_frontfacing", "debug_diffuse", "debug_specular", "debug_roughness"};

auto tracer_funcs = std::vector<ygl::trace_func>{ygl::trace_path,
    ygl::trace_path_volume, ygl::trace_direct, ygl::trace_environment,
    ygl::trace_eyelight, ygl::trace_path_nomis, ygl::trace_path_naive,
    ygl::trace_direct_nomis, ygl::trace_debug_normal, ygl::trace_debug_albedo,
    ygl::trace_debug_texcoord, ygl::trace_debug_frontfacing,
    ygl::trace_debug_diffuse, ygl::trace_debug_specular,
    ygl::trace_debug_roughness};

void draw_widgets(GLFWwindow* win) {
    static auto first_time = true;
    auto app = (app_state*)glfwGetWindowUserPointer(win);
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    if (first_time) {
        ImGui::SetNextWindowPos({0, 0});
        ImGui::SetNextWindowSize({320, 0});
        ImGui::SetNextWindowCollapsed(true);
        first_time = false;
    }
    if (ImGui::Begin("yitrace")) {
        ImGui::LabelText("scene", "%s", app->filename.c_str());
        ImGui::LabelText("image", "%d x %d @ %d", app->img.width,
            app->img.height, app->sample);
        if (ImGui::TreeNode("render settings")) {
            auto cam_names = std::vector<std::string>();
            for (auto cam : app->scn->cameras) cam_names.push_back(cam->name);
            auto edited = 0;
            edited += ImGui::Combo("camera", &app->camid, cam_names);
            edited +=
                ImGui::SliderInt("resolution", &app->resolution, 256, 4096);
            edited += ImGui::SliderInt("nsamples", &app->nsamples, 16, 4096);
            edited += ImGui::Combo("tracer", &app->tracer_id, tracer_names);
            edited += ImGui::SliderInt("nbounces", &app->nbounces, 1, 10);
            edited += ImGui::SliderInt("seed", (int*)&app->seed, 0, 1000);
            edited += ImGui::SliderInt("pratio", &app->pratio, 1, 64);
            if (edited) app->update_list.push_back({"app", app});
            ImGui::LabelText("time/sample", "%0.3lf",
                (app->sample) ? (ygl::get_time() - app->trace_start) /
                                    (1000000000.0 * app->sample) :
                                0.0);
            ImGui::TreePop();
        }
        if (ImGui::TreeNode("view settings")) {
            ImGui::SliderFloat("exposure", &app->exposure, -5, 5);
            ImGui::SliderFloat("gamma", &app->gamma, 1, 3);
            ImGui::ColorEdit4("background", &app->background.x);
            ImGui::SliderFloat("zoom", &app->imscale, 0.1, 10);
            ImGui::Checkbox("zoom to fit", &app->zoom_to_fit);
            ImGui::SameLine();
            ImGui::Checkbox("fps", &app->navigation_fps);
            auto mouse_x = 0.0, mouse_y = 0.0;
            glfwGetCursorPos(win, &mouse_x, &mouse_y);
            auto ij = ygl::get_image_coords(
                ygl::vec2f{(float)mouse_x, (float)mouse_y}, app->imcenter,
                app->imscale, {app->img.width, app->img.height});
            ImGui::DragInt2("mouse", &ij.x);
            if (ij.x >= 0 && ij.x < app->img.width && ij.y >= 0 &&
                ij.y < app->img.height) {
                ImGui::ColorEdit4("pixel", &app->img.at(ij.x, ij.y).x);
            } else {
                auto zero4f_ = ygl::zero4f;
                ImGui::ColorEdit4("pixel", &zero4f_.x);
            }
            ImGui::TreePop();
        }
        if (ImGui::TreeNode("scene tree")) {
            draw_glwidgets_scene_tree(
                "", app->scn, app->selection, app->update_list, 200);
            ImGui::TreePop();
        }
        if (ImGui::TreeNode("scene object")) {
            draw_glwidgets_scene_inspector(
                "", app->scn, app->selection, app->update_list, 200);
            ImGui::TreePop();
        }
    }
    ImGui::End();
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

void draw(GLFWwindow* win) {
    auto app = (app_state*)glfwGetWindowUserPointer(win);
    auto win_width = 0, win_height = 0;
    auto fb_width = 0, fb_height = 0;
    glfwGetFramebufferSize(win, &fb_width, &fb_height);
    glViewport(0, 0, fb_width, fb_height);
    glfwGetWindowSize(win, &win_width, &win_height);
    glClearColor(0.8f, 0.8f, 0.8f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    ygl::center_image4f(app->imcenter, app->imscale,
        {app->display.width, app->display.height}, {win_width, win_height},
        app->zoom_to_fit);
    draw_glimage({app->display.width, app->display.height},
        app->display.pxl.data(), true, {win_width, win_height}, app->imcenter,
        app->imscale);
    draw_widgets(win);
    glfwSwapBuffers(win);
}

bool update(app_state* app) {
    // exit if no updated
    if (app->update_list.empty()) return false;

    // stop renderer
    ygl::trace_async_stop(app->threads, app->stop);

    // update BVH
    for (auto& sel : app->update_list) {
        if (sel.first == "shape") {
            ygl::refit_bvh((ygl::shape*)sel.second);
            ygl::refit_bvh(app->scn);
        }
        if (sel.first == "instance") { ygl::refit_bvh(app->scn); }
        if (sel.first == "node") {
            ygl::update_transforms(app->scn, 0);
            ygl::refit_bvh(app->scn);
        }
    }
    app->update_list.clear();

    auto cam = app->scn->cameras.at(app->camid);
    auto tracer_func = tracer_funcs.at(app->tracer_id);
    app->trace_start = ygl::get_time();
    ygl::trace_async_start(app->scn, cam, app->nsamples, tracer_func, app->img,
        app->display, app->rng, app->threads, app->stop, app->sample,
        app->exposure, app->gamma, app->filmic, app->pratio, app->nbounces,
        app->pixel_clamp, app->seed);

    // updated
    return true;
}

// run ui loop
void run_ui(app_state* app) {
    // window
    auto ww = ygl::clamp(app->img.width, 256, 1440);
    auto wh = ygl::clamp(app->img.height, 256, 1440);
    auto win = init_glfw_window(ww, wh, "yitrace", app, draw);

    // init widgets
    init_imgui_widgets(win);

    // loop
    auto mouse_pos = ygl::zero2f, last_pos = ygl::zero2f;
    while (!glfwWindowShouldClose(win)) {
        last_pos = mouse_pos;
        double mouse_posx, mouse_posy;
        glfwGetCursorPos(win, &mouse_posx, &mouse_posy);
        mouse_pos = {(float)mouse_posx, (float)mouse_posy};
        auto mouse_left =
            glfwGetMouseButton(win, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
        auto mouse_right =
            glfwGetMouseButton(win, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS;
        auto alt_down = glfwGetKey(win, GLFW_KEY_LEFT_ALT) == GLFW_PRESS ||
                        glfwGetKey(win, GLFW_KEY_RIGHT_ALT) == GLFW_PRESS;
        auto shift_down = glfwGetKey(win, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
                          glfwGetKey(win, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS;
        auto widgets_active = ImGui::GetWidgetsActiveExt();

        // handle mouse and keyboard for navigation
        if ((mouse_left || mouse_right) && !alt_down && !widgets_active) {
            auto dolly = 0.0f;
            auto pan = ygl::zero2f;
            auto rotate = ygl::zero2f;
            if (mouse_left && !shift_down)
                rotate = (mouse_pos - last_pos) / 100.0f;
            if (mouse_right) dolly = (mouse_pos.x - last_pos.x) / 100.0f;
            if (mouse_left && shift_down) pan = (mouse_pos - last_pos) / 100.0f;
            auto cam = app->scn->cameras.at(app->camid);
            ygl::camera_turntable(cam->frame, cam->focus, rotate, dolly, pan);
            app->update_list.push_back({"camera", cam});
        }

        // selection
        if ((mouse_left || mouse_right) && alt_down && !widgets_active) {
            auto ij = ygl::get_image_coords(mouse_pos, app->imcenter,
                app->imscale, {app->img.width, app->img.height});
            if (ij.x < 0 || ij.x >= app->img.width || ij.y < 0 ||
                ij.y >= app->img.height) {
                auto cam = app->scn->cameras.at(app->camid);
                auto ray = eval_camera_ray(cam, ij.x, ij.y, app->img.width,
                    app->img.height, {0.5f, 0.5f}, ygl::zero2f);
                auto isec = intersect_ray(app->scn, ray);
                if (isec.ist) app->selection = isec.ist;
            }
        }

        // update
        update(app);

        // draw
        draw(win);

        // event hadling
        if (!(mouse_left || mouse_right) && !widgets_active)
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        glfwPollEvents();
    }
}

int main(int argc, char* argv[]) {
    // application
    auto app = new app_state();

    // parse command line
    auto parser = ygl::make_cmdline_parser(
        argc, argv, "progressive path tracing", "yitrace");
    app->camid = ygl::parse_int(parser, "--camera", 0, "Camera index.");
    app->resolution = ygl::parse_int(
        parser, "--resolution,-r", 512, "Image vertical resolution.");
    app->nsamples =
        ygl::parse_int(parser, "--nsamples,-s", 4096, "Number of samples.");
    app->tracer_id =
        ygl::parse_enum(parser, "--tracer,-t", 0, "Trace type.", tracer_names);
    app->nbounces =
        ygl::parse_int(parser, "--nbounces", 4, "Maximum number of bounces.");
    app->pixel_clamp =
        ygl::parse_int(parser, "--pixel-clamp", 100, "Final pixel clamping.");
    app->seed = ygl::parse_int(
        parser, "--seed", 7, "Seed for the random number generators.");
    auto embree =
        ygl::parse_flag(parser, "--embree", false, "Use Embree ratracer");
    auto double_sided = ygl::parse_flag(
        parser, "--double-sided", false, "Double-sided rendering.");
    auto add_skyenv =
        ygl::parse_flag(parser, "--add-skyenv", false, "Add missing env map");
    auto quiet =
        ygl::parse_flag(parser, "--quiet", false, "Print only errors messages");
    app->imfilename = ygl::parse_string(
        parser, "--output-image,-o", "out.hdr", "Image filename");
    app->filename = ygl::parse_string(
        parser, "scene", "scene.json", "Scene filename", true);
    ygl::check_cmdline(parser);

    // scene loading
    if (!quiet) printf("loading scene %s\n", app->filename.c_str());
    auto load_start = ygl::get_time();
    try {
        app->scn = ygl::load_scene(app->filename);
    } catch (const std::exception& e) {
        printf("cannot load scene %s\n", app->filename.c_str());
        printf("error: %s\n", e.what());
        exit(1);
    }
    if (!quiet)
        printf("loading in %s\n",
            ygl::format_duration(ygl::get_time() - load_start).c_str());

    // tesselate
    if (!quiet) printf("tesselating scene elements\n");
    ygl::tesselate_subdivs(app->scn);

    // add components
    if (!quiet) printf("adding scene elements\n");
    if (add_skyenv && app->scn->environments.empty()) {
        app->scn->environments.push_back(ygl::make_sky_environment("sky"));
        app->scn->textures.push_back(app->scn->environments.back()->ke_txt);
    }
    if (double_sided)
        for (auto mat : app->scn->materials) mat->double_sided = true;
    if (app->scn->cameras.empty())
        app->scn->cameras.push_back(
            ygl::make_bbox_camera("<view>", ygl::compute_bbox(app->scn)));
    ygl::add_missing_names(app->scn);
    for (auto& err : ygl::validate(app->scn))
        printf("warning: %s\n", err.c_str());

    // build bvh
    if (!quiet) printf("building bvh\n");
    auto bvh_start = ygl::get_time();
    ygl::build_bvh(app->scn);
#if YGL_EMBREE
    if (embree) ygl::build_bvh_embree(app->scn);
#endif
    if (!quiet)
        printf("building bvh in %s\n",
            ygl::format_duration(ygl::get_time() - bvh_start).c_str());

    // init renderer
    if (!quiet) printf("initializing lights\n");
    ygl::init_lights(app->scn);

    // fix renderer type if no lights
    if (app->scn->lights.empty() && app->scn->environments.empty() &&
        app->tracer_id != 3) {
        if (!quiet)
            printf("no lights presents, switching to eyelight shader\n");
        app->tracer_id = 3;
    }

    // prepare renderer
    auto cam = app->scn->cameras.at(app->camid);
    auto tracer_func = tracer_funcs.at(app->tracer_id);
    app->img = ygl::make_image4f(ygl::image_width(cam, app->resolution),
        ygl::image_height(cam, app->resolution));
    app->display = ygl::make_image4f(ygl::image_width(cam, app->resolution),
        ygl::image_height(cam, app->resolution));
    app->rng = ygl::make_trace_rngs(ygl::image_width(cam, app->resolution),
        ygl::image_height(cam, app->resolution), app->seed);

    // initialize rendering objects
    if (!quiet) printf("starting async renderer\n");
    app->trace_start = ygl::get_time();
    ygl::trace_async_start(app->scn, cam, app->nsamples, tracer_func, app->img,
        app->display, app->rng, app->threads, app->stop, app->sample,
        app->exposure, app->gamma, app->filmic, app->pratio, app->nbounces,
        app->pixel_clamp, app->seed);

    // run interactive
    run_ui(app);

    // cleanup
    ygl::trace_async_stop(app->threads, app->stop);
    delete app;

    // done
    return 0;
}
