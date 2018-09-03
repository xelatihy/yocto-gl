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
#include "imgui/imgui.h"
#include "imgui/imgui_ext.h"
#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl2.h"

struct app_state {
    std::vector<std::string> filenames;
    std::vector<std::string> names;
    std::vector<ygl::image4f> imgs;
    int img_id = 0;

    // image display
    ygl::image4f display;
    ygl::vec4f pxl_min = ygl::zero4f, pxl_max = ygl::zero4f;
    float lum_min = 0, lum_max = 0;
    float exposure = 0;
    float gamma = 2.2f;
    bool filmic = false;

    // viewing properties
    ygl::vec2f imcenter = ygl::zero2f;
    float imscale = 1;
    bool zoom_to_fit = false;
    ygl::vec4f background = {0.8f, 0.8f, 0.8f, 0};
};

// compute min/max
void update_minmax(app_state* app) {
    app->pxl_min = {ygl::flt_max, ygl::flt_max, ygl::flt_max, ygl::flt_max};
    app->pxl_max = {ygl::flt_min, ygl::flt_min, ygl::flt_min, ygl::flt_min};
    app->lum_min = ygl::flt_max;
    app->lum_max = ygl::flt_min;
    for (auto p : app->imgs.at(app->img_id).pxl) {
        app->pxl_min.x = ygl::min(app->pxl_min.x, p.x);
        app->pxl_min.y = ygl::min(app->pxl_min.y, p.y);
        app->pxl_min.z = ygl::min(app->pxl_min.z, p.z);
        app->pxl_min.w = ygl::min(app->pxl_min.w, p.w);
        app->pxl_max.x = ygl::max(app->pxl_max.x, p.x);
        app->pxl_max.y = ygl::max(app->pxl_max.y, p.y);
        app->pxl_max.z = ygl::max(app->pxl_max.z, p.z);
        app->pxl_max.w = ygl::max(app->pxl_max.w, p.w);
        app->lum_min = ygl::min(app->lum_min, ygl::luminance(xyz(p)));
        app->lum_max = ygl::max(app->lum_max, ygl::luminance(xyz(p)));
    }
}

void update_display_image(app_state* app) {
    app->display = app->imgs.at(app->img_id);
    if (ygl::is_hdr_filename(app->filenames.at(app->img_id))) {
        app->display = ygl::tonemap_image4f(
            app->display, app->exposure, app->gamma, app->filmic);
    }
}

void draw_widgets(GLFWwindow* win) {
    static auto first_time = true;
    auto app = (app_state*)glfwGetWindowUserPointer(win);
    ImGui_ImplOpenGL2_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    if (first_time) {
        ImGui::SetNextWindowPos({0, 0});
        ImGui::SetNextWindowSize({320, 0});
        ImGui::SetNextWindowCollapsed(true);
        first_time = false;
    }
    if (ImGui::Begin("yimview")) {
        auto edited = 0;
        edited += ImGui::Combo("image", &app->img_id, app->names);
        auto& img = app->imgs.at(app->img_id);
        ImGui::LabelText("filename", "%s", app->filenames[app->img_id].c_str());
        ImGui::LabelText("size", "%d x %d ", img.width, img.height);
        edited += ImGui::SliderFloat("exposure", &app->exposure, -5, 5);
        edited += ImGui::SliderFloat("gamma", &app->gamma, 1, 3);
        edited += ImGui::Checkbox("filmic", &app->filmic);
        if (edited) {
            update_display_image(app);
            update_minmax(app);
        }
        ImGui::SliderFloat("zoom", &app->imscale, 0.1, 10);
        ImGui::Checkbox("zoom to fit", &app->zoom_to_fit);
        ImGui::ColorEdit4("background", &app->background.x);
        ImGui::Separator();
        auto mouse_x = 0.0, mouse_y = 0.0;
        glfwGetCursorPos(win, &mouse_x, &mouse_y);
        auto ij =
            ygl::get_image_coords(ygl::vec2f{(float)mouse_x, (float)mouse_y},
                app->imcenter, app->imscale, {img.width, img.height});
        ImGui::DragInt2("mouse", &ij.x);
        auto pixel = ygl::zero4f;
        if (ij.x >= 0 && ij.x < img.width && ij.y >= 0 && ij.y < img.height) {
            pixel = img.at(ij.x, ij.y);
        }
        ImGui::ColorEdit4("pixel", &pixel.x);
        if (!img.pxl.empty()) {
            ImGui::DragFloat4("pxl min", &app->pxl_min.x);
            ImGui::DragFloat4("pxl max", &app->pxl_max.y);
            ImGui::DragFloat("lum min", &app->lum_min);
            ImGui::DragFloat("lum max", &app->lum_max);
        }
    }
    ImGui::End();
    ImGui::Render();
    ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
}

void draw(GLFWwindow* win) {
    auto app = (app_state*)glfwGetWindowUserPointer(win);
    glClearColor(0.8f, 0.8f, 0.8f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glRasterPos2f(-1, 1);
    glPixelZoom(2, -2);
    glDrawPixels(app->display.width, app->display.height, GL_RGBA, GL_FLOAT,
        app->display.pxl.data());
    draw_widgets(win);
    glfwSwapBuffers(win);
}

void run_ui(app_state* app) {
    // window
    auto ww = ygl::min(app->imgs[0].width, 1440);
    auto wh = ygl::min(app->imgs[0].height, 1440);
    if (!glfwInit()) throw std::runtime_error("cannot open glwindow");

    auto win = glfwCreateWindow(ww, wh, "yimview", nullptr, nullptr);
    glfwMakeContextCurrent(win);
    glfwSwapInterval(1);  // Enable vsync

    glfwSetWindowRefreshCallback(win, draw);
    glfwSetWindowUserPointer(win, app);

    // init widgets
    ImGui::CreateContext();
    ImGui::GetIO().IniFilename = nullptr;
    ImGui_ImplGlfw_InitForOpenGL(win, true);
    ImGui_ImplOpenGL2_Init();
    ImGui::StyleColorsDark();

    // window values
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
        auto widgets_active = ImGui::GetWidgetsActiveExt();

        // handle mouse
        if (mouse_left && !widgets_active)
            app->imcenter += mouse_pos - last_pos;
        if (mouse_right && !widgets_active)
            app->imscale *= powf(2, (mouse_pos.x - last_pos.x) * 0.001f);

        // draw
        draw(win);

        // event hadling
        if (mouse_left || mouse_right || widgets_active) {
            glfwPollEvents();
        } else {
            glfwWaitEvents();
        }
    }

    // cleanup
    glfwTerminate();
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
    auto quiet = ygl::parse_flag(
        parser, "--quiet,-q", false, "Print only errors messages");
    app->filenames =
        ygl::parse_strings(parser, "images", {}, "image filenames", true);
    ygl::check_cmdline(parser);

    // loading images
    for (auto filename : app->filenames) {
        if (!quiet) printf("loading %s\n", filename.c_str());
        app->imgs.push_back(ygl::load_image4f(filename));
        app->names.push_back(ygl::get_filename(filename));
    }
    app->img_id = 0;

    // init
    update_display_image(app);
    update_minmax(app);

    // run ui
    run_ui(app);

    // cleanup
    delete app;

    // done
    return 0;
}
