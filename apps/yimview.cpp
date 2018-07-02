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
#include "yglui.h"

#include <memory>
#include <unordered_map>

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
    bool widgets_open = false;
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
        app->pxl_max.x = ygl::min(app->pxl_max.x, p.x);
        app->pxl_max.y = ygl::min(app->pxl_max.y, p.y);
        app->pxl_max.z = ygl::min(app->pxl_max.z, p.z);
        app->pxl_max.w = ygl::min(app->pxl_max.w, p.w);
        app->lum_min = ygl::min(app->lum_min, ygl::luminance(xyz(p)));
        app->lum_max = ygl::max(app->lum_max, ygl::luminance(xyz(p)));
    }
}

void update_display_image(app_state* app) {
    app->display = app->imgs.at(app->img_id);
    if (ygl::is_hdr_filename(app->filenames.at(app->img_id))) {
        app->display = ygl::tonemap_image(
            app->display, app->exposure, app->gamma, app->filmic);
    }
}

void draw_widgets(GLFWwindow* win, app_state* app) {
    if (ygl::begin_widgets_frame(win, "yimview", &app->widgets_open)) {
        auto edited = 0;
        edited += ImGui::Combo("image", &app->img_id, app->names);
        auto& img = app->imgs.at(app->img_id);
        ImGui::LabelText("filename", "%s", app->filenames[app->img_id].c_str());
        ImGui::LabelText(
            "size", "%d x %d ", img.width, img.height);
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
        auto ij = ygl::get_image_coords(
            ygl::vec2f{(float)mouse_x, (float)mouse_y}, app->imcenter,
            app->imscale, {img.width, img.height});
        ImGui::DragInt2("mouse", &ij.x);
        auto pixel = ygl::zero4f;
        if (ij.x >= 0 && ij.x < img.width && ij.y >= 0 &&
            ij.y < img.height) {
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
    ygl::end_widgets_frame();
}

void draw(GLFWwindow* win) {
    auto app = (app_state*)glfwGetWindowUserPointer(win);
    draw_glimage(win, app->display, app->imcenter, app->imscale,
        app->zoom_to_fit, app->background);
    draw_widgets(win, app);
    glfwSwapBuffers(win);
}

void run_ui(app_state* app) {
    // window
    auto ww = ygl::clamp(app->imgs[0].width, 512, 1024);
    auto wh = ygl::clamp(app->imgs[0].height, 512, 1024);
    auto win = ygl::make_window(ww, wh, "yimview", app, draw);

    // init widgets
    ygl::init_widgets(win);

    // window values
    auto mouse_pos = ygl::zero2f, last_pos = ygl::zero2f;
    auto mouse_button = 0;
    while (!glfwWindowShouldClose(win)) {
        last_pos = mouse_pos;
        glfwGetCursorPosExt(win, &mouse_pos.x, &mouse_pos.y);
        mouse_button = glfwGetMouseButtonIndexExt(win);
        auto widgets_active = ImGui::GetWidgetsActiveExt();

        // handle mouse
        if (mouse_button && !widgets_active) {
            if (mouse_button == 1) app->imcenter += mouse_pos - last_pos;
            if (mouse_button == 2) {
                app->imscale *= powf(2, (mouse_pos.x - last_pos.x) * 0.001f);
            }
        }

        // draw
        draw(win);

        // event hadling
        if (mouse_button || widgets_active) {
            glfwPollEvents();
        } else {
            glfwWaitEvents();
        }
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
    auto quiet = ygl::parse_flag(
        parser, "--quiet,-q", false, "Print only errors messages");
    app->filenames =
        ygl::parse_strings(parser, "images", {}, "image filenames", true);
    ygl::check_cmdline(parser);

    // loading images
    for (auto filename : app->filenames) {
        if (!quiet) printf("loading %s\n", filename.c_str());
        app->imgs.push_back(ygl::load_image(filename));
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
