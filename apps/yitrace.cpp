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
    std::shared_ptr<ygl::scene> scn = nullptr;
    std::string filename = "scene.json"s;
    std::string imfilename = "out.obj"s;

    // rendering params
    int camid = 0;                             // camera index
    int resolution = 512;                      // image resolution
    int nsamples = 256;                        // number of samples
    std::string tracer = "pathtrace"s;         // tracer name
    ygl::trace_func tracef = ygl::trace_path;  // tracer
    int nbounces = 8;                          // max depth
    int seed = 7;                              // seed
    float pixel_clamp = 100.0f;                // pixel clamping
    bool double_sided = false;                 // double sided
    bool add_skyenv = false;                   // add sky environment
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
    uint gl_txt = 0;
    uint gl_prog = 0, gl_vbo = 0, gl_ebo;
    ygl::scene_selection selection = {};
    std::vector<ygl::scene_selection> update_list;
    bool quiet = false;
    bool navigation_fps = false;
    bool rendering = false;
};

auto trace_names = std::vector<std::string>{
    "pathtrace",
    "direct",
    "environment",
    "eyelight",
    "pathtrace_nomis",
    "pathtrace_naive",
    "direct_nomis",
    "debug_normal",
    "debug_albedo",
    "debug_diffuse",
    "debug_specular",
    "debug_roughness",
    "debug_texcoord",
    "debug_frontfacing",
};

auto tracer_names = std::unordered_map<std::string, ygl::trace_func>{
    {"pathtrace", ygl::trace_path},
    {"direct", ygl::trace_direct},
    {"environment", ygl::trace_environment},
    {"eyelight", ygl::trace_eyelight},
    {"pathtrace-nomis", ygl::trace_path_nomis},
    {"pathtrace-naive", ygl::trace_path_naive},
    {"direct-nomis", ygl::trace_direct_nomis},
    {"debug-normal", ygl::trace_debug_normal},
    {"debug-albedo", ygl::trace_debug_albedo},
    {"debug-texcoord", ygl::trace_debug_texcoord},
    {"debug-frontfacing", ygl::trace_debug_frontfacing},
};

void draw_image(GLFWwindow* win) {
    auto app = (app_state*)glfwGetWindowUserPointer(win);
    auto window_size = ygl::zero2i, framebuffer_size = ygl::zero2i;
    glfwGetWindowSize(win, &window_size.x, &window_size.y);
    glfwGetFramebufferSize(win, &framebuffer_size.x, &framebuffer_size.y);

    auto& img = app->trace_state.display;
    ygl::center_image(app->imframe, {img.width(), img.height()}, window_size,
        app->zoom_to_fit);

    assert(glGetError() == GL_NO_ERROR);
    glBindTexture(GL_TEXTURE_2D, app->gl_txt);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, img.width(), img.height(), GL_RGBA,
        GL_FLOAT, img.data());
    assert(glGetError() == GL_NO_ERROR);

    assert(glGetError() == GL_NO_ERROR);
    glViewport(0, 0, framebuffer_size.x, framebuffer_size.y);
    glClearColor(app->background.x, app->background.y, app->background.z,
        app->background.w);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_BLEND);
    glUseProgram(app->gl_prog);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, app->gl_txt);
    assert(glGetError() == GL_NO_ERROR);

    assert(glGetError() == GL_NO_ERROR);
    glUniform2f(glGetUniformLocation(app->gl_prog, "win_size"), window_size.x,
        window_size.y);
    glUniform2f(glGetUniformLocation(app->gl_prog, "txt_size"), img.width(),
        img.height());
    glUniformMatrix3x2fv(glGetUniformLocation(app->gl_prog, "frame"), 1, false,
        &app->imframe.x.x);
    glUniform1i(glGetUniformLocation(app->gl_prog, "img"), 0);
    assert(glGetError() == GL_NO_ERROR);

    assert(glGetError() == GL_NO_ERROR);
    glEnableVertexAttribArray(
        glGetAttribLocation(app->gl_prog, "vert_texcoord"));
    glBindBuffer(GL_ARRAY_BUFFER, app->gl_vbo);
    glVertexAttribPointer(glGetAttribLocation(app->gl_prog, "vert_texcoord"), 2,
        GL_FLOAT, false, 0, 0);
    assert(glGetError() == GL_NO_ERROR);

    assert(glGetError() == GL_NO_ERROR);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, app->gl_ebo);
    glDrawElements(GL_TRIANGLES, 3 * 4, GL_UNSIGNED_INT, 0);
    assert(glGetError() == GL_NO_ERROR);

    glUseProgram(0);
    glDisable(GL_BLEND);
}

void draw_widgets(GLFWwindow* win) {
    auto app = (app_state*)glfwGetWindowUserPointer(win);
    if (begin_widgets_frame(win, "yitrace")) {
        ImGui::LabelText("scene", "%s", app->filename.c_str());
        ImGui::LabelText("image", "%d x %d @ %d", app->trace_state.img.width(),
            app->trace_state.img.height(), app->trace_state.sample);
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
                {app->trace_state.img.width(), app->trace_state.img.height()});
            ImGui::DragInt2("mouse", &ij.x);
            if (ij.x >= 0 && ij.x < app->trace_state.img.width() && ij.y >= 0 &&
                ij.y < app->trace_state.img.height()) {
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
    end_widgets_frame();
}

void draw(GLFWwindow* win) {
    draw_image(win);
    draw_widgets(win);
    glfwSwapBuffers(win);
}

static const char* vertex =
    R"(
    #version 330

    layout(location = 0) in vec2 vert_texcoord;

    uniform mat3x2 frame;
    uniform vec2 txt_size;
    uniform vec2 win_size;
    uniform sampler2D img;

    out vec2 texcoord;

    void main() {
        texcoord = vert_texcoord.xy;
        vec2 pos = frame * vec3(txt_size.x * (vert_texcoord.x - 0.5), 
                                txt_size.y * (vert_texcoord.y - 0.5), 1);
        vec2 upos = 2 * pos / win_size - vec2(1,1);
        upos.y = - upos.y;
        gl_Position = vec4(upos.x, upos.y, 0, 1);
    }

    )";

static const char* fragment =
    R"(
    #version 330

    in vec2 texcoord;

    uniform sampler2D img;

    out vec4 color;

    void main() {
        color = texture(img,texcoord);
    }
    )";

void init_drawimage(GLFWwindow* win) {
    auto app = (app_state*)glfwGetWindowUserPointer(win);
    auto& img = app->trace_state.display;
    app->gl_prog = make_glprogram(vertex, fragment);
    auto uv = std::vector<ygl::vec2f>{{0, 0}, {0, 1}, {1, 1}, {1, 0}};
    auto triangles = std::vector<ygl::vec3i>{{0, 1, 2}, {0, 2, 3}};
    assert(glGetError() == GL_NO_ERROR);
    glGenBuffers(1, &app->gl_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, app->gl_vbo);
    glBufferData(GL_ARRAY_BUFFER, 4 * 2 * 4, uv.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glGenBuffers(1, &app->gl_ebo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, app->gl_ebo);
    glBufferData(
        GL_ELEMENT_ARRAY_BUFFER, 2 * 3 * 4, triangles.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    assert(glGetError() == GL_NO_ERROR);
    assert(glGetError() == GL_NO_ERROR);
    glGenTextures(1, &app->gl_txt);
    glBindTexture(GL_TEXTURE_2D, app->gl_txt);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img.width(), img.height(), 0,
        GL_RGBA, GL_FLOAT, img.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glBindTexture(GL_TEXTURE_2D, 0);
    assert(glGetError() == GL_NO_ERROR);
}

bool update(const std::shared_ptr<app_state>& app) {
    // exit if no updated
    if (app->update_list.empty()) return false;

    // stop renderer
    ygl::trace_async_stop(app->trace_state);

    // update BVH
    for (auto sel : app->update_list) {
        if (sel.as<ygl::shape>()) {
            if (!app->quiet) std::cout << "refit shape bvh\n";
            ygl::refit_bvh(sel.as<ygl::shape>());
        }
        if (sel.as<ygl::instance>()) {
            if (!app->quiet) std::cout << "refit scene bvh\n";
            ygl::refit_bvh(app->scn, false);
        }
        if (sel.as<ygl::node>()) {
            if (!app->quiet) std::cout << "refit scene bvh\n";
            ygl::update_transforms(app->scn, 0);
            ygl::refit_bvh(app->scn, false);
        }
    }
    app->update_list.clear();

    app->trace_state = {};
    ygl::trace_async_start(app->trace_state, app->scn, app->camid,
        app->resolution, app->nsamples, app->tracef, app->exposure, app->gamma,
        app->filmic, app->pratio, app->nbounces, app->pixel_clamp, app->seed);

    // updated
    return true;
}

// run ui loop
void run_ui(const std::shared_ptr<app_state>& app) {
    // window
    auto win_width = ygl::clamp(app->trace_state.img.width(), 512, 1024);
    auto win_height = ygl::clamp(app->trace_state.img.height(), 512, 1024);
    auto win = make_window(win_width, win_height, "yitrace", app.get(), draw);

    // load textures
    init_drawimage(win);

    // init widget
    init_widgets(win);

    // loop
    auto mouse_pos = ygl::zero2f, last_pos = ygl::zero2f;
    auto mouse_center = ygl::zero3f;
    auto mouse_button = 0, last_button = 0;
    while (!glfwWindowShouldClose(win)) {
        last_pos = mouse_pos;
        last_button = mouse_button;
        glfwGetCursorPosExt(win, &mouse_pos.x, &mouse_pos.y);
        mouse_button = glfwGetMouseButtonIndexExt(win);
        auto alt_down = glfwGetAltKeyExt(win);
        auto shift_down = glfwGetShiftKeyExt(win);
        auto widgets_active = ImGui::GetWidgetsActiveExt();

        // handle mouse and keyboard for navigation
        if (mouse_button && !alt_down && !widgets_active) {
            if (!last_button) {
                if (app->selection.as<ygl::instance>()) {
                    auto ist = app->selection.as<ygl::instance>();
                    mouse_center = transform_point(ist->frame,
                        (ist->shp->bbox.min + ist->shp->bbox.max) / 2);
                } else {
                    mouse_center =
                        (app->scn->bbox.min + app->scn->bbox.max) / 2;
                }
            }
            auto dolly = 0.0f;
            auto pan = ygl::zero2f;
            auto rotate = ygl::zero2f;
            if (mouse_button == 1) rotate = (mouse_pos - last_pos) / 100.0f;
            if (mouse_button == 2) dolly = (mouse_pos.x - last_pos.x) / 100.0f;
            if (mouse_button == 3 || (mouse_button == 1 && shift_down))
                pan = (mouse_pos - last_pos) / 100.0f;
            auto cam = app->scn->cameras.at(app->camid);
            ygl::camera_turntable(
                cam->frame.o, mouse_center, cam->frame.z, rotate, dolly, pan);
            cam->frame = ygl::lookat_frame(
                cam->frame.o, mouse_center, ygl::vec3f{0, 1, 0});
            app->update_list.push_back(cam);
        }

        // selection
        if (mouse_button && alt_down && !widgets_active) {
            auto ij = ygl::get_image_coords(mouse_pos, app->imframe,
                ygl::vec2i{app->trace_state.img.width(),
                    app->trace_state.img.height()});
            if (ij.x < 0 || ij.x >= app->trace_state.img.width() || ij.y < 0 ||
                ij.y >= app->trace_state.img.height()) {
                auto cam = app->scn->cameras.at(app->camid);
                auto ray = eval_camera_ray(cam, ij,
                    ygl::vec2i{app->trace_state.img.width(),
                        app->trace_state.img.height()},
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
        glfwPollEvents();
    }
}

int main(int argc, char* argv[]) {
    // create empty scene
    auto app = std::make_shared<app_state>();

    // parse command line
    CLI::App parser("progressive path tracing", "yitrace");
    parser.add_option("--camera", app->camid, "Camera index.");
    parser.add_option(
        "--resolution,-r", app->resolution, "Image vertical resolution.");
    parser.add_option("--nsamples,-s", app->nsamples, "Number of samples.");
    parser.add_option("--tracer,-t", app->tracer, "Trace type.")
        ->check([](const std::string& s) -> std::string {
            if (tracer_names.find(s) == tracer_names.end())
                throw CLI::ValidationError("unknown tracer name");
            return s;
        });
    parser.add_option(
        "--nbounces", app->nbounces, "Maximum number of bounces.");
    parser.add_option(
        "--pixel-clamp", app->pixel_clamp, "Final pixel clamping.");
    parser.add_option(
        "--seed", app->seed, "Seed for the random number generators.");
    parser.add_option("--exposure,-e", app->exposure, "Hdr exposure");
    parser.add_flag(
        "--double-sided,-D", app->double_sided, "Double-sided rendering.");
    parser.add_flag("--add-skyenv,-E", app->add_skyenv, "add missing env map");
    parser.add_flag("--quiet,-q", app->quiet, "Print only errors messages");
    parser.add_option(
        "--pration", app->pratio, "Preview ratio for async rendering.");
    parser.add_option("--output-image,-o", app->imfilename, "Image filename");
    parser.add_option("scene", app->filename, "Scene filename")->required(true);
    try {
        parser.parse(argc, argv);
    } catch (const CLI::ParseError& e) { return parser.exit(e); }
    app->tracef = tracer_names.at(app->tracer);

    // scene loading
    if (!app->quiet) std::cout << "loading scene" << app->filename << "\n";
    auto load_start = ygl::get_time();
    try {
        app->scn = ygl::load_scene(app->filename);
    } catch (const std::exception& e) {
        std::cout << "cannot load scene " << app->filename << "\n";
        std::cout << "error: " << e.what() << "\n";
        exit(1);
    }
    if (!app->quiet)
        std::cout << "loading in "
                  << ygl::format_duration(ygl::get_time() - load_start) << "\n";

    // tesselate
    if (!app->quiet) std::cout << "tesselating scene elements\n";
    ygl::update_tesselation(app->scn);

    // update bbox and transforms
    ygl::update_transforms(app->scn);
    ygl::update_bbox(app->scn);

    // add components
    if (!app->quiet) std::cout << "adding scene elements\n";
    if (app->add_skyenv && app->scn->environments.empty()) {
        app->scn->environments.push_back(ygl::make_sky_environment("sky"));
        app->scn->textures.push_back(app->scn->environments.back()->ke_txt);
    }
    if (app->double_sided)
        for (auto mat : app->scn->materials) mat->double_sided = true;
    if (app->scn->cameras.empty())
        app->scn->cameras.push_back(
            ygl::make_bbox_camera("<view>", app->scn->bbox));
    ygl::add_missing_names(app->scn);
    for (auto err : ygl::validate(app->scn))
        std::cout << "warning: " << err << "\n";

    // build bvh
    if (!app->quiet) std::cout << "building bvh\n";
    auto bvh_start = ygl::get_time();
    ygl::update_bvh(app->scn);
    if (!app->quiet)
        std::cout << "building bvh in "
                  << ygl::format_duration(ygl::get_time() - bvh_start) << "\n";

    // init renderer
    if (!app->quiet) std::cout << "initializing lights\n";
    ygl::update_lights(app->scn);

    // fix renderer type if no lights
    if (app->scn->lights.empty() && app->scn->environments.empty() &&
        app->tracer != "eyelight") {
        if (!app->quiet)
            std::cout << "no lights presents, switching to eyelight shader\n";
        app->tracer = "eyelight";
        app->tracef = ygl::trace_eyelight;
    }

    // initialize rendering objects
    if (!app->quiet) std::cout << "starting async renderer\n";
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
