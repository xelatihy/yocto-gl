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
#ifdef _WIN32
#include <windows.h>
#endif
#include <GL/gl.h>
#endif
#include <GLFW/glfw3.h>

struct app_state {
    // image
    std::string name = "";
    std::string filename = "";
    ygl::image4f img = {};

    // image adjustment
    float exposure = 0;
    float gamma = 2.2f;
    bool filmic = false;
    bool is_hdr = false;

    // display image
    ygl::image4f display;
};

void update_display_image(app_state* app) {
    app->display = app->img;
    if (app->is_hdr) {
        app->display = ygl::tonemap_image4f(
            app->display, app->exposure, app->gamma, app->filmic);
    }
}

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
    update_display_image(app);
}

void run_ui(app_state* app) {
    // window
    auto ww = ygl::clamp(app->img.width, 256, 1440);
    auto wh = ygl::clamp(app->img.height, 256, 1440);
    if (!glfwInit()) throw std::runtime_error("cannot open glwindow");

    auto win = glfwCreateWindow(ww, wh, "yimview", nullptr, nullptr);
    glfwMakeContextCurrent(win);
    glfwSwapInterval(1);  // Enable vsync

    glfwSetCharCallback(win, char_callback);
    glfwSetWindowRefreshCallback(win, draw);
    glfwSetWindowUserPointer(win, app);

    // window values
    while (!glfwWindowShouldClose(win)) {
        draw(win);
        glfwWaitEvents();
    }
}

int main(int argc, char* argv[]) {
    // prepare application
    auto app = new app_state();

    // command line params
    auto parser =
        ygl::make_cmdline_parser(argc, argv, "view images", "yimview");
    app->gamma = ygl::parse_float(parser, "--gamma,-g", 2.2f, "display gamma");
    app->exposure =
        ygl::parse_float(parser, "--exposure,-e", 0, "display exposure");
    app->filmic = ygl::parse_flag(parser, "--filmic", false, "display filmic");
    app->filename = ygl::parse_string(
        parser, "image", app->filename, "image filename", true);
    ygl::check_cmdline(parser);

    // loading image
    app->name = ygl::get_filename(app->filename);
    app->is_hdr = ygl::is_hdr_filename(app->filename);
    app->img = ygl::load_image4f(app->filename);

    // prepare display image
    update_display_image(app);

    // run ui
    run_ui(app);

    // cleanup
    delete app;

    // done
    return 0;
}
