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

struct image_stats {
    ygl::bbox4f pxl_bounds = {ygl::zero4f, ygl::zero4f};
    ygl::bbox1f lum_bounds = {0, 0};
};

struct app_image {
    // original data
    std::string filename;
    std::string name;
    ygl::image4f img;
    bool is_hdr = false;

    // diplay image
    ygl::image4f display;
    ygl::uint gl_txt = 0;

    // image stats
    image_stats stats;

    // tonemapping values
    float exposure = 0;
    float gamma = 2.2f;
    bool filmic = false;

    // computation futures
    std::atomic<bool> load_done, display_done, stats_done, texture_done;
    std::thread load_thread, display_thread, stats_thread;
    std::atomic<bool> display_stop;
    std::string error_msg = "";

    // viewing properties
    ygl::vec2f imcenter = ygl::zero2f;
    float imscale = 1;
    bool zoom_to_fit = true;

    // cleanup and stop threads
    ~app_image() {
        display_stop = true;
        if (load_thread.joinable()) load_thread.join();
        if (display_thread.joinable()) display_thread.join();
        if (stats_thread.joinable()) stats_thread.join();
    }
};

struct app_state {
    // images
    std::vector<app_image*> imgs;
    int img_id = 0;

    ~app_state() {
        for (auto img : imgs) delete img;
    }
};

// compute min/max
void update_image_stats(app_image* img) {
    img->stats_done = false;
    img->stats.pxl_bounds = ygl::invalid_bbox4f;
    img->stats.lum_bounds = ygl::invalid_bbox1f;
    for (auto p : img->img.pxl) {
        img->stats.pxl_bounds += p;
        img->stats.lum_bounds += ygl::luminance(xyz(p));
    }
    img->stats_done = true;
}

void update_display_image(app_image* img) {
    auto start = ygl::get_time();
    img->display_done = false;
    img->texture_done = false;
    if (img->is_hdr) {
        if (img->img.width * img->img.height > 1024 * 1024) {
            auto nthreads = std::thread::hardware_concurrency();
            auto threads = std::vector<std::thread>();
            for (auto tid = 0; tid < nthreads; tid++) {
                threads.push_back(std::thread([img, tid, nthreads]() {
                    for (auto j = tid; j < img->img.height; j += nthreads) {
                        if (img->display_stop) break;
                        for (auto i = 0; i < img->img.width; i++) {
                            img->display.at(i, j) =
                                ygl::tonemap_exposuregamma(img->img.at(i, j),
                                    img->exposure, img->gamma, img->filmic);
                        }
                    }
                }));
            }
            for (auto& t : threads) t.join();
        } else {
            for (auto j = 0; j < img->img.height; j++) {
                if (img->display_stop) break;
                for (auto i = 0; i < img->img.width; i++) {
                    img->display.at(i, j) =
                        ygl::tonemap_exposuregamma(img->img.at(i, j),
                            img->exposure, img->gamma, img->filmic);
                }
            }
        }
    } else {
        img->display = img->img;
    }
    img->display_done = true;
    if (!img->display_stop) {
        printf("display: %s\n",
            ygl::format_duration(ygl::get_time() - start).c_str());
        fflush(stdout);
    }
}

// load image
void load_image(app_image* img) {
    auto start = ygl::get_time();
    try {
        img->img = ygl::load_image4f(img->filename);
    } catch (...) {
        img->error_msg = "cannot load image";
        return;
    }
    img->load_done = true;
    img->display = img->img;
    printf("load: %s\n", ygl::format_duration(ygl::get_time() - start).c_str());
    fflush(stdout);
    img->display_thread = std::thread(update_display_image, img);
    img->stats_thread = std::thread(update_image_stats, img);
}

void draw_glwidgets(ygl::glwindow* win) {
    auto app = (app_state*)ygl::get_user_pointer(win);
    ygl::begin_glwidgets_frame(win);
    if (ygl::begin_glwidgets_window(win, "yimview")) {
        auto edited = 0;
        edited +=
            ygl::draw_combobox_glwidget(win, "image", app->img_id, app->imgs);
        auto img = app->imgs.at(app->img_id);
        ygl::draw_imgui_label(win, "filename", "%s", img->filename.c_str());
        auto status = std::string();
        if (img->error_msg != "")
            status = "error: " + img->error_msg;
        else if (!img->load_done)
            status = "loading...";
        else if (!img->display_done)
            status = "displaying...";
        else
            status = "done";
        ygl::draw_imgui_label(win, "status", status.c_str());
        ygl::draw_imgui_label(
            win, "size", "%d x %d ", img->img.width, img->img.height);
        edited +=
            ygl::draw_slider_glwidget(win, "exposure", img->exposure, -5, 5);
        edited += ygl::draw_slider_glwidget(win, "gamma", img->gamma, 1, 3);
        edited += ygl::draw_checkbox_glwidget(win, "filmic", img->filmic);
        if (edited) {
            if (img->display_thread.joinable()) {
                img->display_stop = true;
                img->display_thread.join();
            }
            img->display_stop = false;
            img->display_thread = std::thread(update_display_image, img);
        }
        ygl::draw_slider_glwidget(win, "zoom", img->imscale, 0.1, 10);
        ygl::draw_checkbox_glwidget(win, "zoom to fit", img->zoom_to_fit);
        ygl::draw_separator_glwidget(win);
        auto mouse_pos = ygl::get_glmouse_pos(win);
        auto ij = ygl::get_image_coords(mouse_pos, img->imcenter, img->imscale,
            {img->img.width, img->img.height});
        ygl::draw_dragger_glwidget(win, "mouse", ij);
        auto pixel = ygl::zero4f;
        if (ij.x >= 0 && ij.x < img->img.width && ij.y >= 0 &&
            ij.y < img->img.height) {
            pixel = img->img.at(ij.x, ij.y);
        }
        ygl::draw_coloredit_glwidget(win, "pixel", pixel);
        auto stats = (img->stats_done) ? img->stats : image_stats{};
        ygl::draw_dragger_glwidget(win, "pxl min", stats.pxl_bounds.min);
        ygl::draw_dragger_glwidget(win, "pxl max", stats.pxl_bounds.max);
        ygl::draw_dragger_glwidget(win, "lum min", stats.lum_bounds.min);
        ygl::draw_dragger_glwidget(win, "lum max", stats.lum_bounds.max);
    }
    ygl::end_glwidgets_frame(win);
}

void draw(ygl::glwindow* win) {
    auto app = (app_state*)ygl::get_user_pointer(win);
    auto img = app->imgs.at(app->img_id);
    auto win_size = ygl::get_glwindow_size(win);
    auto fb_size = ygl::get_glframebuffer_size(win);
    ygl::set_glviewport(fb_size);
    ygl::clear_glframebuffer(ygl::vec4f{0.8f, 0.8f, 0.8f, 1.0f});
    if (img->gl_txt) {
        ygl::center_image4f(img->imcenter, img->imscale,
            {img->display.width, img->display.height}, win_size,
            img->zoom_to_fit);
        ygl::draw_glimage(img->gl_txt,
            {img->display.width, img->display.height}, win_size, img->imcenter,
            img->imscale);
    }
    draw_glwidgets(win);
    ygl::swap_glbuffers(win);
}

void update(app_state* app) {
    for (auto img : app->imgs) {
        if (!img->display_done || img->texture_done) continue;
        if (!img->gl_txt) {
            img->gl_txt = ygl::make_gltexture(img->display, false, false);
        } else {
            ygl::update_gltexture(img->gl_txt, img->display, false, false);
        }
        img->texture_done = true;
    }
}

void run_ui(app_state* app) {
    // window
    auto img = app->imgs.at(app->img_id);
    auto ww = ygl::clamp(img->img.width, 512, 1440);
    auto wh = ygl::clamp(img->img.height, 512, 1440);
    auto win = ygl::make_glwindow(ww, wh, "yimview", app, draw);

    // init widgets
    ygl::init_glwidgets(win);

    // center image
    ygl::center_image4f(img->imcenter, img->imscale,
        {img->img.width, img->img.height}, {ww, wh}, false);

    // window values
    auto mouse_pos = ygl::zero2f, last_pos = ygl::zero2f;
    while (!ygl::should_glwindow_close(win)) {
        last_pos = mouse_pos;
        mouse_pos = ygl::get_glmouse_pos(win);
        auto mouse_left = ygl::get_glmouse_left(win);
        auto mouse_right = ygl::get_glmouse_right(win);
        auto widgets_active = ygl::get_glwidgets_active(win);

        // handle mouse
        if (mouse_left && !widgets_active)
            img->imcenter += mouse_pos - last_pos;
        if (mouse_right && !widgets_active)
            img->imscale *= powf(2, (mouse_pos.x - last_pos.x) * 0.001f);

        // update
        update(app);

        // draw
        draw(win);

        // event hadling
        // could also wait if (mouse_left || mouse_right || widgets_active)
        ygl::process_glevents(win);
    }

    // cleanup
    ygl::delete_glwindow(win);
}

int main(int argc, char* argv[]) {
    // prepare application
    auto app = new app_state();

    // command line params
    auto parser =
        ygl::make_cmdline_parser(argc, argv, "view images", "yimview");
    auto gamma = ygl::parse_float(parser, "--gamma,-g", 2.2f, "display gamma");
    auto exposure =
        ygl::parse_float(parser, "--exposure,-e", 0, "display exposure");
    auto filmic = ygl::parse_flag(parser, "--filmic", false, "display filmic");
    // auto quiet = ygl::parse_flag(
    //     parser, "--quiet,-q", false, "Print only errors messages");
    auto filenames =
        ygl::parse_strings(parser, "images", {}, "image filenames", true);
    ygl::check_cmdline(parser);

    // loading images
    for (auto filename : filenames) {
        auto img = new app_image();
        img->filename = filename;
        img->name = ygl::get_filename(filename);
        img->is_hdr = ygl::is_hdr_filename(filename);
        img->exposure = exposure;
        img->gamma = gamma;
        img->filmic = filmic;
        img->load_thread = std::thread(load_image, img);
        app->imgs.push_back(img);
    }
    app->img_id = 0;

    // wait a bit to see if the first image gets loaded
    std::this_thread::sleep_for(std::chrono::milliseconds(100));

    // run ui
    run_ui(app);

    // cleanup
    delete app;

    // done
    return 0;
}
