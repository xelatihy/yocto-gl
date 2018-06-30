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
    ygl::vec4f pxl_min = ygl::zero4f, pxl_max = ygl::zero4f;
    float lum_min = 0, lum_max = 0;

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
    std::vector<gimage*> imgs;
    gimage* img = nullptr;

    // viewing properties
    uint gl_prog = 0, gl_pbo = 0, gl_tbo = 0, gl_ebo = 0;
    bool widgets_open = false;
    ygl::vec2f imcenter = ygl::zero2f;
    float imscale = 1;
    bool zoom_to_fit = false;
    ygl::vec4f background = {0.8f, 0.8f, 0.8f, 0};

    ~app_state() {
        for (auto v : imgs) delete v;
    }
};

// compute min/max
void update_minmax(gimage* img) {
    img->pxl_min = {ygl::flt_max, ygl::flt_max, ygl::flt_max, ygl::flt_max};
    img->pxl_max = {ygl::flt_min, ygl::flt_min, ygl::flt_min, ygl::flt_min};
    img->lum_min = ygl::flt_max;
    img->lum_max = ygl::flt_min;
    for (auto p : img->img.pxl) {
        img->pxl_min.x = ygl::min(img->pxl_min.x, p.x);
        img->pxl_min.y = ygl::min(img->pxl_min.y, p.y);
        img->pxl_min.z = ygl::min(img->pxl_min.z, p.z);
        img->pxl_min.w = ygl::min(img->pxl_min.w, p.w);
        img->pxl_max.x = ygl::min(img->pxl_max.x, p.x);
        img->pxl_max.y = ygl::min(img->pxl_max.y, p.y);
        img->pxl_max.z = ygl::min(img->pxl_max.z, p.z);
        img->pxl_max.w = ygl::min(img->pxl_max.w, p.w);
        img->lum_min = ygl::min(img->lum_min, ygl::luminance(xyz(p)));
        img->lum_max = ygl::max(img->lum_max, ygl::luminance(xyz(p)));
    }
}

// Loads a generic image
gimage* load_gimage(
    const std::string& filename, float exposure, float gamma, bool filmic) {
    auto img = new gimage();
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

gimage* diff_gimage(const gimage* a, const gimage* b, bool color) {
    if (a->img.width != b->img.width || a->img.height != b->img.height)
        return nullptr;
    auto d = new gimage();
    d->name = "diff " + a->name + " " + b->name;
    d->filename = "";
    d->img = a->img;
    if (color) {
        for (auto j = 0; j < a->img.height; j++) {
            for (auto i = 0; i < a->img.width; i++) {
                d->img.at(i, j) = {
                    std::abs(a->img.at(i, j).x - b->img.at(i, j).x),
                    std::abs(a->img.at(i, j).y - b->img.at(i, j).y),
                    std::abs(a->img.at(i, j).z - b->img.at(i, j).z),
                    std::max(a->img.at(i, j).w, b->img.at(i, j).w)};
            }
        }
    } else {
        for (auto j = 0; j < a->img.height; j++) {
            for (auto i = 0; i < a->img.width; i++) {
                auto la = (a->img.at(i, j).x + a->img.at(i, j).y +
                              a->img.at(i, j).z) /
                          3;
                auto lb = (b->img.at(i, j).x + b->img.at(i, j).y +
                              b->img.at(i, j).z) /
                          3;
                auto ld = fabsf(la - lb);
                d->img.at(i, j) = {
                    ld, ld, ld, std::max(a->img.at(i, j).w, b->img.at(i, j).w)};
            }
        }
    }
    update_minmax(d);
    return d;
}

void update_display_image(gimage* img) {
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
            "size", "%d x %d ", app->img->img.width, app->img->img.height);
        auto edited = 0;
        edited += ImGui::SliderFloat("exposure", &app->img->exposure, -5, 5);
        edited += ImGui::SliderFloat("gamma", &app->img->gamma, 1, 3);
        edited += ImGui::Checkbox("filmic", &app->img->filmic);
        if (edited) app->img->updated = true;
        ImGui::SliderFloat("zoom", &app->imscale, 0.1, 10);
        ImGui::Checkbox("zoom to fit", &app->zoom_to_fit);
        ImGui::ColorEdit4("background", &app->background.x);
        ImGui::Separator();
        auto mouse_x = 0.0, mouse_y = 0.0;
        glfwGetCursorPos(win, &mouse_x, &mouse_y);
        auto ij = ygl::get_image_coords(
            ygl::vec2f{(float)mouse_x, (float)mouse_y}, app->imcenter,
            app->imscale, {app->img->img.width, app->img->img.height});
        ImGui::DragInt2("mouse", &ij.x);
        auto pixel = ygl::zero4f;
        if (ij.x >= 0 && ij.x < app->img->img.width && ij.y >= 0 &&
            ij.y < app->img->img.height) {
            pixel = app->img->img.at(ij.x, ij.y);
        }
        ImGui::ColorEdit4("pixel", &pixel.x);
        if (!app->img->img.pxl.empty()) {
            ImGui::DragFloat4("pxl min", &app->img->pxl_min.x);
            ImGui::DragFloat4("pxl max", &app->img->pxl_max.y);
            ImGui::DragFloat("lum min", &app->img->lum_min);
            ImGui::DragFloat("lum max", &app->img->lum_max);
        }
    }
    ygl::end_widgets_frame();
}

void draw(GLFWwindow* win) {
    auto app = (app_state*)glfwGetWindowUserPointer(win);
    draw_glimage(win, app->img->display, app->imcenter, app->imscale,
        app->zoom_to_fit, app->background);
    draw_widgets(win, app);
    glfwSwapBuffers(win);
}

void run_ui(const std::shared_ptr<app_state>& app) {
    // window
    auto ww = ygl::clamp(app->imgs[0]->img.width, 512, 1024);
    auto wh = ygl::clamp(app->imgs[0]->img.height, 512, 1024);
    auto win = ygl::make_window(ww, wh, "yimview", app.get(), draw);

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
