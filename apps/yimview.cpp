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
using namespace ygl;

struct image_stats {
    bbox4f pxl_bounds = {zero4f, zero4f};
    bbox1f lum_bounds = {0, 0};
};

struct app_image {
    // original data
    std::string filename;
    std::string name;
    image<vec4f> img;
    bool is_hdr = false;

    // diplay image
    image<vec4f> display;
    uint gl_txt = 0;

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
    vec2f imcenter = zero2f;
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
    img->stats.pxl_bounds = invalid_bbox4f;
    img->stats.lum_bounds = invalid_bbox1f;
    for (auto p : img->img) {
        img->stats.pxl_bounds += p;
        img->stats.lum_bounds += luminance(xyz(p));
    }
    img->stats_done = true;
}

void update_display_image(app_image* img) {
    auto start = get_time();
    img->display_done = false;
    img->texture_done = false;
    if (img->is_hdr) {
        if (img->img.size().x * img->img.size().y > 1024 * 1024) {
            auto nthreads = std::thread::hardware_concurrency();
            auto threads = std::vector<std::thread>();
            for (auto tid = 0; tid < nthreads; tid++) {
                threads.push_back(std::thread([img, tid, nthreads]() {
                    for (auto j = tid; j < img->img.size().y; j += nthreads) {
                        if (img->display_stop) break;
                        for (auto i = 0; i < img->img.size().x; i++) {
                            img->display[{i, j}] =
                                tonemap_exposuregamma(img->img[{i, j}],
                                    img->exposure, img->gamma, img->filmic);
                        }
                    }
                }));
            }
            for (auto& t : threads) t.join();
        } else {
            for (auto j = 0; j < img->img.size().y; j++) {
                if (img->display_stop) break;
                for (auto i = 0; i < img->img.size().x; i++) {
                    img->display[{i, j}] =
                        tonemap_exposuregamma(img->img[{i, j}], img->exposure,
                            img->gamma, img->filmic);
                }
            }
        }
    } else {
        img->display = img->img;
    }
    img->display_done = true;
    if (!img->display_stop) {
        printf("display: %s\n", format_duration(get_time() - start).c_str());
        fflush(stdout);
    }
}

// load image
void load_image(app_image* img) {
    auto start = get_time();
    try {
        img->img = load_image4f(img->filename);
    } catch (...) {
        img->error_msg = "cannot load image";
        return;
    }
    img->load_done = true;
    img->display = img->img;
    printf("load: %s\n", format_duration(get_time() - start).c_str());
    fflush(stdout);
    img->display_thread = std::thread(update_display_image, img);
    img->stats_thread = std::thread(update_image_stats, img);
}

void draw_glwidgets(glwindow* win) {
    auto app = (app_state*)get_user_pointer(win);
    begin_glwidgets_frame(win);
    if (begin_glwidgets_window(win, "yimview")) {
        auto edited = 0;
        edited += draw_combobox_glwidget(win, "image", app->img_id, app->imgs);
        auto img = app->imgs.at(app->img_id);
        draw_imgui_label(win, "filename", "%s", img->filename.c_str());
        auto status = std::string();
        if (img->error_msg != "")
            status = "error: " + img->error_msg;
        else if (!img->load_done)
            status = "loading...";
        else if (!img->display_done)
            status = "displaying...";
        else
            status = "done";
        draw_imgui_label(win, "status", status.c_str());
        draw_imgui_label(
            win, "size", "%d x %d ", img->img.size().x, img->img.size().y);
        edited += draw_slider_glwidget(win, "exposure", img->exposure, -5, 5);
        edited += draw_slider_glwidget(win, "gamma", img->gamma, 1, 3);
        edited += draw_checkbox_glwidget(win, "filmic", img->filmic);
        if (edited) {
            if (img->display_thread.joinable()) {
                img->display_stop = true;
                img->display_thread.join();
            }
            img->display_stop = false;
            img->display_thread = std::thread(update_display_image, img);
        }
        draw_slider_glwidget(win, "zoom", img->imscale, 0.1, 10);
        draw_checkbox_glwidget(win, "zoom to fit", img->zoom_to_fit);
        draw_separator_glwidget(win);
        auto mouse_pos = get_glmouse_pos(win);
        auto ij = get_image_coords(mouse_pos, img->imcenter, img->imscale,
            {img->img.size().x, img->img.size().y});
        draw_dragger_glwidget(win, "mouse", ij);
        auto pixel = zero4f;
        if (ij.x >= 0 && ij.x < img->img.size().x && ij.y >= 0 &&
            ij.y < img->img.size().y) {
            pixel = img->img[ij];
        }
        draw_coloredit_glwidget(win, "pixel", pixel);
        auto stats = (img->stats_done) ? img->stats : image_stats{};
        draw_dragger_glwidget(win, "pxl min", stats.pxl_bounds.min);
        draw_dragger_glwidget(win, "pxl max", stats.pxl_bounds.max);
        draw_dragger_glwidget(win, "lum min", stats.lum_bounds.min);
        draw_dragger_glwidget(win, "lum max", stats.lum_bounds.max);
    }
    end_glwidgets_frame(win);
}

void draw(glwindow* win) {
    auto app = (app_state*)get_user_pointer(win);
    auto img = app->imgs.at(app->img_id);
    auto win_size = get_glwindow_size(win);
    auto fb_size = get_glframebuffer_size(win);
    set_glviewport(fb_size);
    clear_glframebuffer(vec4f{0.8f, 0.8f, 0.8f, 1.0f});
    if (img->gl_txt) {
        center_image4f(img->imcenter, img->imscale,
            {img->display.size().x, img->display.size().y}, win_size,
            img->zoom_to_fit);
        draw_glimage(img->gl_txt,
            {img->display.size().x, img->display.size().y}, win_size,
            img->imcenter, img->imscale);
    }
    draw_glwidgets(win);
    swap_glbuffers(win);
}

void update(app_state* app) {
    for (auto img : app->imgs) {
        if (!img->display_done || img->texture_done) continue;
        if (!img->gl_txt) {
            img->gl_txt = make_gltexture(img->display, false, false);
        } else {
            update_gltexture(img->gl_txt, img->display, false, false);
        }
        img->texture_done = true;
    }
}

void run_ui(app_state* app) {
    // window
    auto img = app->imgs.at(app->img_id);
    auto win_size = clamp(img->img.size(), 512, 1440);
    auto win = make_glwindow(win_size, "yimview", app, draw);

    // init widgets
    init_glwidgets(win);

    // center image
    center_image4f(
        img->imcenter, img->imscale, img->img.size(), win_size, false);

    // window values
    auto mouse_pos = zero2f, last_pos = zero2f;
    while (!should_glwindow_close(win)) {
        last_pos = mouse_pos;
        mouse_pos = get_glmouse_pos(win);
        auto mouse_left = get_glmouse_left(win);
        auto mouse_right = get_glmouse_right(win);
        auto widgets_active = get_glwidgets_active(win);

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
        process_glevents(win);
    }

    // cleanup
    delete_glwindow(win);
}

int main(int argc, char* argv[]) {
    // prepare application
    auto app = new app_state();

    // command line params
    auto parser = make_cmdline_parser(argc, argv, "view images", "yimview");
    auto gamma = parse_arg(parser, "--gamma,-g", 2.2f, "display gamma");
    auto exposure = parse_arg(parser, "--exposure,-e", 0.0f, "display exposure");
    auto filmic = parse_arg(parser, "--filmic", false, "display filmic");
    // auto quiet = parse_flag(
    //     parser, "--quiet,-q", false, "Print only errors messages");
    auto filenames =
        parse_args(parser, "images", std::vector<std::string>{}, "image filenames", true);
    check_cmdline(parser);

    // loading images
    for (auto filename : filenames) {
        auto img = new app_image();
        img->filename = filename;
        img->name = get_filename(filename);
        img->is_hdr = is_hdr_filename(filename);
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
