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

#include <memory>
#include <unordered_map>

// Generic image that contains either an HDR or an LDR image, giving access
// to both. This is helpful when writing viewers or generic image
// manipulation code
struct gimage {
    std::string name = ""s;
    std::string filename = ""s;  // path
    ygl::image4f img = {};       // image

    // min/max values
    ygl::bbox4f pxl_bounds = ygl::invalid_bbox4f;
    ygl::bbox1f lum_bounds = ygl::invalid_bbox1f;

    // image adjustment
    bool updated = true;
    float exposure = 0;
    float gamma = 2.2f;
    bool filmic = false;
    bool is_hdr = false;

    // display image
    ygl::image4f display;
};

struct app_state {
    std::vector<std::shared_ptr<gimage>> imgs;
    std::shared_ptr<gimage> img = nullptr;

    // viewing properties
    uint gl_prog = 0, gl_pbo = 0, gl_tbo = 0, gl_ebo = 0;
    bool widgets_open = false;
    ygl::frame2f imframe = ygl::identity_frame2f;
    bool zoom_to_fit = false;
    ygl::vec4f background = {0.8f, 0.8f, 0.8f, 0};
};

// compute min/max
void update_minmax(std::shared_ptr<gimage> img) {
    img->pxl_bounds = ygl::invalid_bbox4f;
    img->lum_bounds = ygl::invalid_bbox1f;
    for (auto& p : img->img.pxl) {
        img->pxl_bounds += p;
        img->lum_bounds += ygl::luminance(xyz(p));
    }
}

// Loads a generic image
std::shared_ptr<gimage> load_gimage(
    const std::string& filename, float exposure, float gamma, bool filmic) {
    auto img = std::make_shared<gimage>();
    img->filename = filename;
    img->name = ygl::get_filename(filename);
    img->exposure = exposure;
    img->gamma = gamma;
    img->filmic = filmic;
    img->is_hdr = ygl::is_hdr_filename(filename);
    img->img = ygl::load_image(filename);
    update_minmax(img);
    return img;
}

std::shared_ptr<gimage> diff_gimage(
    std::shared_ptr<gimage> a, std::shared_ptr<gimage> b, bool color) {
    if (a->img.size != b->img.size) return nullptr;
    auto d = std::make_shared<gimage>();
    d->name = "diff " + a->name + " " + b->name;
    d->filename = "";
    d->img = a->img;
    if (color) {
        for (auto j = 0; j < a->img.size.y; j++) {
            for (auto i = 0; i < a->img.size.x; i++) {
                d->img[{i, j}] = {std::abs(a->img[{i, j}].x - b->img[{i, j}].x),
                    std::abs(a->img[{i, j}].y - b->img[{i, j}].y),
                    std::abs(a->img[{i, j}].z - b->img[{i, j}].z),
                    std::max(a->img[{i, j}].w, b->img[{i, j}].w)};
            }
        }
    } else {
        for (auto j = 0; j < a->img.size.y; j++) {
            for (auto i = 0; i < a->img.size.x; i++) {
                auto la =
                    (a->img[{i, j}].x + a->img[{i, j}].y + a->img[{i, j}].z) /
                    3;
                auto lb =
                    (b->img[{i, j}].x + b->img[{i, j}].y + b->img[{i, j}].z) /
                    3;
                auto ld = fabsf(la - lb);
                d->img[{i, j}] = {
                    ld, ld, ld, std::max(a->img[{i, j}].w, b->img[{i, j}].w)};
            }
        }
    }
    update_minmax(d);
    return d;
}

void update_display_image(const std::shared_ptr<gimage>& img) {
    img->display = img->img;
    if (img->is_hdr) {
        img->display = ygl::tonemap_image(
            img->display, img->exposure, img->gamma, img->filmic);
    }
}

void draw_widgets(GLFWwindow* win, app_state* app) {
    if (ygl::begin_widgets_frame(win, "yimview", &app->widgets_open)) {
        ImGui::Combo("image", &app->img, app->imgs, false);
        ImGui::LabelText("filename", "%s", app->img->filename.c_str());
        ImGui::LabelText(
            "size", "%d x %d ", app->img->img.size.x, app->img->img.size.y);
        auto edited = 0;
        edited += ImGui::SliderFloat("exposure", &app->img->exposure, -5, 5);
        edited += ImGui::SliderFloat("gamma", &app->img->gamma, 1, 3);
        edited += ImGui::Checkbox("filmic", &app->img->filmic);
        if (edited) app->img->updated = true;
        auto zoom = app->imframe.x.x;
        if (ImGui::SliderFloat("zoom", &zoom, 0.1, 10))
            app->imframe.x.x = app->imframe.y.y = zoom;
        ImGui::Checkbox("zoom to fit", &app->zoom_to_fit);
        ImGui::ColorEdit4("background", &app->background.x);
        ImGui::Separator();
        auto mouse_x = 0.0, mouse_y = 0.0;
        glfwGetCursorPos(win, &mouse_x, &mouse_y);
        auto ij =
            ygl::get_image_coords(ygl::vec2f{(float)mouse_x, (float)mouse_y},
                app->imframe, {app->img->img.size.x, app->img->img.size.y});
        ImGui::DragInt2("mouse", &ij.x);
        auto pixel = ygl::zero4f;
        if (ij.x >= 0 && ij.x < app->img->img.size.x && ij.y >= 0 &&
            ij.y < app->img->img.size.y) {
            pixel = app->img->img.at(ij);
        }
        ImGui::ColorEdit4("pixel", &pixel.x);
        if (!app->img->img.pxl.empty()) {
            ImGui::DragFloat4("pxl min", &app->img->pxl_bounds.min.x);
            ImGui::DragFloat4("pxl max", &app->img->pxl_bounds.max.y);
            ImGui::DragFloat("lum min", &app->img->lum_bounds.min);
            ImGui::DragFloat("lum max", &app->img->lum_bounds.max);
        }
    }
    ygl::end_widgets_frame();
}

void draw(GLFWwindow* win) {
    auto app = (app_state*)glfwGetWindowUserPointer(win);
    draw_glimage(win, app->img->display, app->imframe, app->zoom_to_fit,
        app->background);
    draw_widgets(win, app);
    glfwSwapBuffers(win);
}

void run_ui(const std::shared_ptr<app_state>& app) {
    // window
    auto win_size = ygl::clamp(app->imgs[0]->img.size, 512, 1024);
    auto win = make_window(win_size, "yimview", app.get(), draw);

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
            if (mouse_button == 1) app->imframe.o += mouse_pos - last_pos;
            if (mouse_button == 2) {
                auto zoom = powf(2, (mouse_pos.x - last_pos.x) * 0.001f);
                app->imframe =
                    app->imframe * ygl::frame2f{{zoom, 0}, {0, zoom}, {0, 0}};
            }
        }

        // update texture
        if (app->img->updated) {
            update_display_image(app->img);
            app->img->updated = false;
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
    // command line parameters
    std::vector<std::string> filenames = {"img.png"};
    float gamma = 2.2f;     // hdr gamma
    float exposure = 0.0f;  // hdr exposure
    bool filmic = false;    // filmic tone mapping
    bool diff = false;      // compute diffs
    bool lum_diff = false;  // compute luminance diffs
    bool quiet = false;     // quiet mode

    // command line params
    CLI::App parser("view images", "yimview");
    parser.add_option("--gamma,-g", gamma, "display gamma");
    parser.add_option("--exposure,-e", exposure, "display exposure");
    parser.add_flag("--filmic", filmic, "display filmic");
    parser.add_flag("--diff,-d", diff, "compute diff images");
    parser.add_flag("--luminance-diff,-D", lum_diff, "compute luminance diffs");
    parser.add_flag("--quiet,-q", quiet, "Print only errors messages");
    parser.add_option("images", filenames, "image filenames")->required(true);
    try {
        parser.parse(argc, argv);
    } catch (const CLI::ParseError& e) { return parser.exit(e); }

    // prepare application
    auto app = std::make_shared<app_state>();

    // loading images
    for (auto filename : filenames) {
        if (!quiet) std::cout << "loading " << filename << "\n";
        app->imgs.push_back(load_gimage(filename, exposure, gamma, filmic));
    }
    app->img = app->imgs.at(0);
    if (diff) {
        for (auto i = 0; i < app->imgs.size(); i++) {
            for (auto j = i + 1; j < app->imgs.size(); j++) {
                if (!quiet)
                    std::cout << "diffing " << app->imgs[i]->filename << " "
                              << app->imgs[j]->filename << "\n";
                auto diff = diff_gimage(app->imgs[i], app->imgs[j], !lum_diff);
                if (diff) app->imgs.push_back(diff);
            }
        }
    }

    // run ui
    run_ui(app);

    // done
    return 0;
}
