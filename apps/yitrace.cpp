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
    ygl::bvh_tree* bvh = nullptr;

    // rendering params
    std::string filename = "scene.json";
    std::string imfilename = "out.obj";
    ygl::trace_params prm = {};

    // rendering state
    ygl::trace_state* stt = nullptr;

    // view image
    ygl::vec2f imcenter = ygl::zero2f;
    float imscale = 1;
    bool zoom_to_fit = true;
    ygl::vec4f background = {0.8f, 0.8f, 0.8f, 0};
    bool widgets_open = false;
    void* selection = nullptr;
    std::vector<std::pair<std::string, void*>> update_list;
    bool navigation_fps = false;
    bool quiet = false;
    int64_t trace_start = 0;

    ~app_state() {
        if (scn) delete scn;
        if (bvh) delete bvh;
        if (stt) delete stt;
    }
};

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
        ImGui::LabelText("image", "%d x %d @ %d", app->stt->img.width,
            app->stt->img.height, app->stt->sample);
        if (ImGui::TreeNode("render settings")) {
            auto cam_names = std::vector<std::string>();
            for (auto cam : app->scn->cameras) cam_names.push_back(cam->name);
            auto edited = 0;
            edited += ImGui::Combo("camera", &app->prm.camid, cam_names);
            edited += ImGui::SliderInt(
                "resolution", &app->prm.yresolution, 256, 4096);
            edited +=
                ImGui::SliderInt("nsamples", &app->prm.nsamples, 16, 4096);
            edited += ImGui::Combo(
                "tracer", (int*)&app->prm.tracer, ygl::trace_type_names);
            edited += ImGui::SliderInt("nbounces", &app->prm.nbounces, 1, 10);
            edited += ImGui::SliderInt("seed", (int*)&app->prm.seed, 0, 1000);
            edited +=
                ImGui::SliderInt("pratio", &app->prm.preview_ratio, 1, 64);
            if (edited) app->update_list.push_back({"app", app});
            ImGui::LabelText("time/sample", "%0.3lf",
                (app->stt->sample) ? (ygl::get_time() - app->trace_start) /
                                         (1000000000.0 * app->stt->sample) :
                                     0.0);
            ImGui::TreePop();
        }
        if (ImGui::TreeNode("view settings")) {
            ImGui::SliderFloat("exposure", &app->prm.exposure, -5, 5);
            ImGui::SliderFloat("gamma", &app->prm.gamma, 1, 3);
            ImGui::ColorEdit4("background", &app->background.x);
            ImGui::SliderFloat("zoom", &app->imscale, 0.1, 10);
            ImGui::Checkbox("zoom to fit", &app->zoom_to_fit);
            ImGui::SameLine();
            ImGui::Checkbox("fps", &app->navigation_fps);
            auto mouse_x = 0.0, mouse_y = 0.0;
            glfwGetCursorPos(win, &mouse_x, &mouse_y);
            auto ij = ygl::get_image_coords(
                ygl::vec2f{(float)mouse_x, (float)mouse_y}, app->imcenter,
                app->imscale, {app->stt->img.width, app->stt->img.height});
            ImGui::DragInt2("mouse", &ij.x);
            if (ij.x >= 0 && ij.x < app->stt->img.width && ij.y >= 0 &&
                ij.y < app->stt->img.height) {
                ImGui::ColorEdit4("pixel", &app->stt->img.at(ij.x, ij.y).x);
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
        {app->stt->display.width, app->stt->display.height},
        {win_width, win_height}, app->zoom_to_fit);
    draw_glimage({app->stt->display.width, app->stt->display.height},
        app->stt->display.pxl.data(), true, {win_width, win_height},
        app->imcenter, app->imscale);
    draw_widgets(win);
    glfwSwapBuffers(win);
}

bool update(app_state* app) {
    // exit if no updated
    if (app->update_list.empty()) return false;

    // stop renderer
    ygl::trace_async_stop(app->stt);

    // update BVH
    for (auto& sel : app->update_list) {
        if (sel.first == "shape") {
            for(auto sid = 0; sid < app->scn->shapes.size(); sid++) {
                if(app->scn->shapes[sid] == sel.second) {
                    ygl::refit_bvh((ygl::shape*)sel.second, app->bvh->shape_bvhs[sid]);
                    break;
                }
            }
            ygl::refit_bvh(app->scn, app->bvh);
        }
        if (sel.first == "instance") { ygl::refit_bvh(app->scn, app->bvh); }
        if (sel.first == "node") {
            ygl::update_transforms(app->scn, 0);
            ygl::refit_bvh(app->scn, app->bvh);
        }
    }
    app->update_list.clear();

    delete app->stt;
    app->trace_start = ygl::get_time();
    app->stt = ygl::make_trace_state(app->scn, app->prm);
    ygl::trace_async_start(app->stt, app->scn, app->bvh, app->prm);

    // updated
    return true;
}

// run ui loop
void run_ui(app_state* app) {
    // window
    auto ww = ygl::clamp(app->stt->img.width, 256, 1440);
    auto wh = ygl::clamp(app->stt->img.height, 256, 1440);
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
            auto cam = app->scn->cameras.at(app->prm.camid);
            ygl::camera_turntable(cam->frame, cam->focus, rotate, dolly, pan);
            app->update_list.push_back({"camera", cam});
        }

        // selection
        if ((mouse_left || mouse_right) && alt_down && !widgets_active) {
            auto ij = ygl::get_image_coords(mouse_pos, app->imcenter,
                app->imscale, {app->stt->img.width, app->stt->img.height});
            if (ij.x < 0 || ij.x >= app->stt->img.width || ij.y < 0 ||
                ij.y >= app->stt->img.height) {
                auto cam = app->scn->cameras.at(app->prm.camid);
                auto ray = eval_camera_ray(cam, ij.x, ij.y, app->stt->img.width,
                    app->stt->img.height, {0.5f, 0.5f}, ygl::zero2f);
                auto isec = intersect_ray(app->scn, app->bvh, ray);
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
    app->prm.camid = ygl::parse_int(parser, "--camera", 0, "Camera index.");
    app->prm.yresolution = ygl::parse_int(
        parser, "--resolution,-r", 512, "Image vertical resolution.");
    app->prm.nsamples =
        ygl::parse_int(parser, "--nsamples,-s", 4096, "Number of samples.");
    app->prm.tracer = (ygl::trace_type)ygl::parse_enum(
        parser, "--tracer,-t", 0, "Tracer type.", ygl::trace_type_names);
    app->prm.nbounces =
        ygl::parse_int(parser, "--nbounces", 4, "Maximum number of bounces.");
    app->prm.pixel_clamp =
        ygl::parse_int(parser, "--pixel-clamp", 100, "Final pixel clamping.");
    app->prm.seed = ygl::parse_int(
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
    app->bvh = ygl::build_bvh(app->scn, true, embree);
    if (!quiet)
        printf("building bvh in %s\n",
            ygl::format_duration(ygl::get_time() - bvh_start).c_str());

    // init renderer
    if (!quiet) printf("initializing lights\n");
    ygl::init_lights(app->scn);

    // fix renderer type if no lights
    if (app->scn->lights.empty() && app->scn->environments.empty() &&
        app->prm.tracer != ygl::trace_type::eyelight) {
        if (!quiet)
            printf("no lights presents, switching to eyelight shader\n");
        app->prm.tracer = ygl::trace_type::eyelight;
    }

    // prepare renderer
    app->stt = ygl::make_trace_state(app->scn, app->prm);

    // initialize rendering objects
    if (!quiet) printf("starting async renderer\n");
    app->trace_start = ygl::get_time();
    ygl::trace_async_start(app->stt, app->scn, app->bvh, app->prm);

    // run interactive
    run_ui(app);

    // cleanup
    ygl::trace_async_stop(app->stt);
    delete app;

    // done
    return 0;
}
