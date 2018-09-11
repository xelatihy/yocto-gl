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
#endif
#include <GLFW/glfw3.h>

// Application state
struct app_state {
    // scene
    std::string filename = "scene.json";
    std::string imfilename = "out.obj";
    ygl::scene* scn = nullptr;
    ygl::trace_params prm = {};
    ygl::trace_state* stt = nullptr;
    bool embree = false;  // Embree acceleration

    // view image
    float exposure = 0;
    float gamma = 2.2f;
    bool filmic = false;

    ~app_state() {
        if (scn) delete scn;
        if (stt) delete stt;
    }
};

void draw(GLFWwindow* win) {
    auto app = (app_state*)glfwGetWindowUserPointer(win);
    glClearColor(0.8f, 0.8f, 0.8f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glRasterPos2f(-1, 1);
    glPixelZoom(2, -2);
    glDrawPixels(app->stt->display.width, app->stt->display.height, GL_RGBA,
        GL_FLOAT, app->stt->display.pxl.data());
    glfwSwapBuffers(win);
}

void restart(app_state* app) {
    // stop renderer
    ygl::trace_async_stop(app->stt);
    delete app->stt;
    app->stt = ygl::make_trace_state(app->scn, app->prm);
    ygl::trace_async_start(app->stt, app->scn, app->prm);
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
            app->prm.camid = (app->prm.camid + 1) % app->scn->cameras.size();
            break;
        case 'C':
            app->prm.camid = (app->prm.camid - 1 + app->scn->cameras.size()) %
                             app->scn->cameras.size();
            break;
        case 't':
            app->prm.tracer = (ygl::trace_type)(
                ((int)app->prm.tracer + 1) % ygl::trace_type_names.size());
            break;
        case 'T':
            app->prm.tracer = (ygl::trace_type)(
                ((int)app->prm.tracer - 1 + ygl::trace_type_names.size()) %
                ygl::trace_type_names.size());
            break;
        case ' ': break;  // restart
    }
    auto cam = app->scn->cameras.at(app->prm.camid);
    printf("trace with %s from %s camera\n",
        ygl::trace_type_names.at((int)app->prm.tracer).c_str(),
        cam->name.c_str());
    printf("tonemap with %g exposure, %g gamma and %d filmic\n", app->exposure,
        app->gamma, (int)app->filmic);
    restart(app);
}

// run ui loop
void run_ui(app_state* app) {
    // window
    auto ww = ygl::clamp(app->stt->img.width, 256, 1440);
    auto wh = ygl::clamp(app->stt->img.height, 256, 1440);
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
            auto cam = app->scn->cameras.at(app->prm.camid);
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
    app->prm.camid = ygl::parse_int(parser, "--camera", 0, "Camera index.");
    app->prm.yresolution = ygl::parse_int(
        parser, "--resolution,-r", 512, "Image vertical resolution.");
    app->prm.nsamples =
        ygl::parse_int(parser, "--nsamples,-s", 4096, "Number of samples.");
    app->prm.tracer = (ygl::trace_type)ygl::parse_enum(
        parser, "--tracer,-t", 0, "Trace type.", ygl::trace_type_names);
    app->prm.nbounces =
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
    app->stt = ygl::make_trace_state(app->scn, app->prm);

    // initialize rendering objects
    printf("starting async renderer\n");
    ygl::trace_async_start(app->stt, app->scn, app->prm);

    // run interactive
    run_ui(app);

    // cleanup
    ygl::trace_async_stop(app->stt);
    delete app;

    // done
    return 0;
}
