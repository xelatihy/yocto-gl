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
using namespace std::literals;

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/glew.h>
#endif
#include <GLFW/glfw3.h>

// Application state
struct app_state {
    // scene
    ygl::scene* scn = nullptr;
    ygl::camera* cam = nullptr;

    // rendering params
    std::string filename = "scene.json"s;
    std::string imfilename = "out.obj"s;
    int resolution = 512;                      // image resolution
    int nsamples = 256;                        // number of samples
    std::string tracer = "pathtrace"s;         // tracer name
    ygl::trace_func tracef = ygl::trace_path;  // tracer
    int nbounces = 4;                          // max depth
    int seed = ygl::trace_default_seed;        // seed
    float pixel_clamp = 100.0f;                // pixel clamping
    int pratio = 8;                            // preview ratio

    // rendering state
    ygl::image4f img = {};
    ygl::image4f display = {};
    std::vector<ygl::rng_state> rng = {};
    bool stop = false;
    std::vector<std::thread> threads;
    int sample = 0;

    // view image
    float exposure = 0;
    float gamma = 2.2f;
    bool filmic = false;
    bool navigation_fps = false;

    ~app_state() {
        if (scn) delete scn;
    }
};

auto trace_names = std::vector<std::string>{"pathtrace", "direct",
    "environment", "eyelight", "pathtrace-nomis", "pathtrace-naive",
    "direct-nomis", "debug_normal", "debug_albedo", "debug_texcoord",
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

void draw(GLFWwindow* win) {
    auto app = (app_state*)glfwGetWindowUserPointer(win);
    glClearColor(0.8f, 0.8f, 0.8f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glRasterPos2f(-1, 1);
    glPixelZoom(2, -2);
    glDrawPixels(app->display.width, app->display.height, GL_RGBA, GL_FLOAT,
        app->display.pxl.data());
    glfwSwapBuffers(win);
}

void restart(app_state* app) {
    // stop renderer
    ygl::trace_async_stop(&app->threads, &app->stop);
    app->tracef = tracer_names.at(app->tracer);
    ygl::trace_async_start(app->scn, app->cam, app->nsamples, app->tracef,
        &app->img, &app->display, &app->rng, &app->threads, &app->stop,
        &app->sample, &app->exposure, &app->gamma, &app->filmic, app->pratio,
        app->nbounces, app->pixel_clamp, app->seed);
}

void char_callback(GLFWwindow* win, unsigned int key) {
    auto app = (app_state*)glfwGetWindowUserPointer(win);
    switch ((char)key) {
        case '1':
            app->exposure = 0;
            app->gamma = 1;
            break;
        case '2':
            app->exposure = 0;
            app->gamma = 2.2f;
            break;
        case '[': app->exposure -= 1; break;
        case ']': app->exposure += 1; break;
        case '{': app->gamma -= 0.1f; break;
        case '}': app->gamma += 0.1f; break;
        case 'f': app->filmic = !app->filmic; break;
    }
    printf("tonemap with %g exposure, %g gamma and %d filmic\n", app->exposure,
        app->gamma, (int)app->filmic);
    restart(app);
}

// run ui loop
void run_ui(app_state* app) {
    // window
    auto ww = ygl::clamp(app->img.width, 512, 1024);
    auto wh = ygl::clamp(app->img.height, 512, 1024);
    if (!glfwInit()) throw std::runtime_error("cannot open glwindow");

    auto win = glfwCreateWindow(ww, wh, "yimview", nullptr, nullptr);
    glfwMakeContextCurrent(win);
    glfwSwapInterval(1);  // Enable vsync

    glfwSetCharCallback(win, char_callback);
    glfwSetWindowRefreshCallback(win, draw);
    glfwSetWindowUserPointer(win, app);

    // loop
    auto mouse_pos = ygl::zero2f, last_pos = ygl::zero2f;
    auto mouse_button = 0;
    while (!glfwWindowShouldClose(win)) {
        auto mx = 0.0, my = 0.0;
        last_pos = mouse_pos;
        glfwGetCursorPos(win, &mx, &my);
        mouse_pos = {(float)mx, (float)my};
        if (glfwGetMouseButton(win, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS)
            mouse_button = 1;
        else if (glfwGetMouseButton(win, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS)
            mouse_button = 2;
        else if (glfwGetMouseButton(win, GLFW_MOUSE_BUTTON_MIDDLE) ==
                 GLFW_PRESS)
            mouse_button = 3;
        else
            mouse_button = 0;
        auto shift_down = glfwGetKey(win, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
                          glfwGetKey(win, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS;

        // handle mouse and keyboard for navigation
        if (mouse_button) {
            auto dolly = 0.0f;
            auto pan = ygl::zero2f;
            auto rotate = ygl::zero2f;
            if (mouse_button == 1) rotate = (mouse_pos - last_pos) / 100.0f;
            if (mouse_button == 2) dolly = (mouse_pos.x - last_pos.x) / 100.0f;
            if (mouse_button == 3 || (mouse_button == 1 && shift_down))
                pan = (mouse_pos - last_pos) / 100.0f;
            ygl::camera_turntable(
                app->cam->frame, app->cam->focus, rotate, dolly, pan);
            restart(app);
        }

        // draw
        draw(win);

        // event hadling
        if (!mouse_button)
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
    auto embree = false;            // use embree
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
    parser.add_flag("--embree", embree, "Use Embree ratracer");
    parser.add_option("--output-image,-o", imfilename, "Image filename");
    parser.add_option("scene", filename, "Scene filename")->required(true);
    try {
        parser.parse(argc, argv);
    } catch (const CLI::ParseError& e) { return parser.exit(e); }

    // scene loading
    auto scn = (ygl::scene*)nullptr;
    if (!quiet) printf("loading scene %s\n", filename.c_str());
    auto load_start = ygl::get_time();
    try {
        scn = ygl::load_scene(filename);
    } catch (const std::exception& e) {
        printf("cannot load scene %s\n", filename.c_str());
        printf("error: %s\n", e.what());
        exit(1);
    }
    if (!quiet)
        printf("loading in %s\n",
            ygl::format_duration(ygl::get_time() - load_start).c_str());

    // tesselate
    if (!quiet) printf("tesselating scene elements\n");
    ygl::tesselate_subdivs(scn);

    // add components
    if (!quiet) printf("adding scene elements\n");
    if (add_skyenv && scn->environments.empty()) {
        scn->environments.push_back(ygl::make_sky_environment("sky"));
        scn->textures.push_back(scn->environments.back()->ke_txt);
    }
    if (double_sided)
        for (auto mat : scn->materials) mat->double_sided = true;
    if (scn->cameras.empty())
        scn->cameras.push_back(
            ygl::make_bbox_camera("<view>", ygl::compute_bbox(scn)));
    ygl::add_missing_names(scn);
    for (auto& err : ygl::validate(scn)) printf("warning: %s\n", err.c_str());

    // build bvh
    if (!quiet) printf("building bvh\n");
    auto bvh_start = ygl::get_time();
    ygl::build_bvh(scn);
#if YGL_EMBREE
    if (embree) ygl::build_bvh_embree(scn);
#endif
    if (!quiet)
        printf("building bvh in %s\n",
            ygl::format_duration(ygl::get_time() - bvh_start).c_str());

    // init renderer
    if (!quiet) printf("initializing lights\n");
    ygl::init_lights(scn);

    // fix renderer type if no lights
    if (scn->lights.empty() && scn->environments.empty() &&
        tracer != "eyelight") {
        if (!quiet)
            printf("no lights presents, switching to eyelight shader\n");
        tracer = "eyelight";
    }

    // prepare application
    auto app = new app_state();
    app->scn = scn;
    app->filename = filename;
    app->imfilename = imfilename;
    app->cam = scn->cameras.at(camid);
    app->resolution = resolution;
    app->nsamples = nsamples;
    app->tracer = tracer;
    app->tracef = tracer_names.at(tracer);
    app->img = ygl::make_image4f(ygl::image_width(app->cam, resolution),
        ygl::image_height(app->cam, resolution));
    app->display = ygl::make_image4f(ygl::image_width(app->cam, resolution),
        ygl::image_height(app->cam, resolution));
    app->rng = ygl::make_trace_rngs(ygl::image_width(app->cam, resolution),
        ygl::image_height(app->cam, resolution), seed);
    app->nbounces = nbounces;
    app->seed = seed;
    app->pixel_clamp = pixel_clamp;
    app->pratio = pratio;

    // initialize rendering objects
    if (!quiet) printf("starting async renderer\n");
    ygl::trace_async_start(app->scn, app->cam, app->nsamples, app->tracef,
        &app->img, &app->display, &app->rng, &app->threads, &app->stop,
        &app->sample, &app->exposure, &app->gamma, &app->filmic, app->pratio,
        app->nbounces, app->pixel_clamp, app->seed);

    // run interactive
    run_ui(app);

    // cleanup
    ygl::trace_async_stop(&app->threads, &app->stop);
    delete app;

    // done
    return 0;
}
