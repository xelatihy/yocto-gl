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

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif
#include <GLFW/glfw3.h>

// Application state
struct app_state {
    // scene
    std::string filename = "scene.json";
    std::string imfilename = "out.obj";
    ygl::scene* scn = nullptr;

    // rendering params
    int camid = 0;         // camera id
    int resolution = 512;  // image resolution
    int nsamples = 256;    // number of samples
    int tracer_id = 0;     // tracer id
    int nbounces = 4;      // max depth
    bool embree = false;   // Embree acceleration

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

    ~app_state() {
        if (scn) delete scn;
    }
};

auto tracer_names = std::vector<std::string>{"pathtrace", "direct",
    "environment", "eyelight", "pathtrace-nomis", "pathtrace-naive",
    "direct-nomis", "debug_normal", "debug_albedo", "debug_texcoord",
    "debug_frontfacing", "debug_diffuse", "debug_specular", "debug_roughness"};

auto tracer_funcs = std::vector<ygl::trace_func>{ygl::trace_path,
    ygl::trace_direct, ygl::trace_environment, ygl::trace_eyelight,
    ygl::trace_path_nomis, ygl::trace_path_naive, ygl::trace_direct_nomis,
    ygl::trace_debug_normal, ygl::trace_debug_albedo, ygl::trace_debug_texcoord,
    ygl::trace_debug_frontfacing, ygl::trace_debug_diffuse,
    ygl::trace_debug_specular, ygl::trace_debug_roughness};

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
    ygl::trace_async_stop(app->threads, app->stop);
    auto cam = app->scn->cameras.at(app->camid);
    auto tracer_func = tracer_funcs.at(app->tracer_id);
    ygl::trace_async_start(app->scn, cam, app->nsamples, tracer_func, app->img,
        app->display, app->rng, app->threads, app->stop, app->sample,
        app->exposure, app->gamma, app->filmic);
}

void char_callback(GLFWwindow* win, unsigned int key) {
    auto app = (app_state*)glfwGetWindowUserPointer(win);
    switch ((char)key) {
        case '[': app->exposure -= 1; break;
        case ']': app->exposure += 1; break;
        case '{': app->gamma -= 0.1f; break;
        case '}': app->gamma += 0.1f; break;
        case 'f': app->filmic = !app->filmic; break;
        case 'c':
            app->camid = (app->camid + 1) % app->scn->cameras.size();
            break;
        case 'C':
            app->camid = (app->camid - 1 + app->scn->cameras.size()) %
                         app->scn->cameras.size();
            break;
        case 't':
            app->tracer_id = (app->tracer_id + 1) % tracer_names.size();
            break;
        case 'T':
            app->tracer_id = (app->tracer_id - 1 + tracer_names.size()) %
                             tracer_names.size();
            break;
        case ' ': break;  // restart
    }
    auto cam = app->scn->cameras.at(app->camid);
    printf("trace with %s from %s camera\n",
        tracer_names.at(app->tracer_id).c_str(), cam->name.c_str());
    printf("tonemap with %g exposure, %g gamma and %d filmic\n", app->exposure,
        app->gamma, (int)app->filmic);
    restart(app);
}

// run ui loop
void run_ui(app_state* app) {
    // window
    auto ww = ygl::clamp(app->img.width, 960, 1440);
    auto wh = ygl::clamp(app->img.height, 960, 1440);
    if (!glfwInit()) throw std::runtime_error("cannot open glwindow");

    auto win = glfwCreateWindow(ww, wh, "yimview", nullptr, nullptr);
    glfwMakeContextCurrent(win);
    glfwSwapInterval(1);  // Enable vsync

    glfwSetCharCallback(win, char_callback);
    glfwSetWindowRefreshCallback(win, draw);
    glfwSetWindowUserPointer(win, app);

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
        auto shift_down = glfwGetKey(win, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
                          glfwGetKey(win, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS;

        // handle mouse and keyboard for navigation
        if (mouse_left || mouse_right) {
            auto dolly = 0.0f;
            auto pan = ygl::zero2f;
            auto rotate = ygl::zero2f;
            if (mouse_left && !shift_down)
                rotate = (mouse_pos - last_pos) / 100.0f;
            if (mouse_right) dolly = (mouse_pos.x - last_pos.x) / 100.0f;
            if (mouse_left && shift_down) pan = (mouse_pos - last_pos) / 100.0f;
            auto cam = app->scn->cameras.at(app->camid);
            ygl::camera_turntable(cam->frame, cam->focus, rotate, dolly, pan);
            restart(app);
        }

        // draw
        draw(win);

        // event hadling
        if (!(mouse_left || mouse_right))
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
    app->embree =
        ygl::parse_flag(parser, "--embree", false, "Use Embree ratracer");
    app->imfilename = ygl::parse_string(
        parser, "--output-image,-o", "out.hdr", "Image filename");
    app->filename = ygl::parse_string(
        parser, "scene", "scene.json", "Scene filename", true);
    ygl::check_cmdline(parser);

    // scene loading
    printf("loading scene %s\n", app->filename.c_str());
    app->scn = ygl::load_scene(app->filename);

    // tesselate
    printf("tesselating scene elements\n");
    ygl::tesselate_subdivs(app->scn);

    // build bvh
    printf("building bvh\n");
    ygl::build_bvh(app->scn);
#if YGL_EMBREE
    if (app->embree) ygl::build_bvh_embree(app->scn);
#endif

    // init renderer
    printf("initializing lights\n");
    ygl::init_lights(app->scn);

    // prepare application
    auto cam = app->scn->cameras.at(app->camid);
    auto tracer_func = tracer_funcs.at(app->tracer_id);
    app->img = ygl::make_image4f(ygl::image_width(cam, app->resolution),
        ygl::image_height(cam, app->resolution));
    app->display = ygl::make_image4f(ygl::image_width(cam, app->resolution),
        ygl::image_height(cam, app->resolution));
    app->rng = ygl::make_trace_rngs(ygl::image_width(cam, app->resolution),
        ygl::image_height(cam, app->resolution));

    // initialize rendering objects
    printf("starting async renderer\n");
    ygl::trace_async_start(app->scn, cam, app->nsamples, tracer_func, app->img,
        app->display, app->rng, app->threads, app->stop, app->sample,
        app->exposure, app->gamma, app->filmic);

    // run interactive
    run_ui(app);

    // cleanup
    ygl::trace_async_stop(app->threads, app->stop);
    delete app;

    // done
    return 0;
}
