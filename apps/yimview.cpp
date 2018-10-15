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
    string       filename;  // filename
    string       outname;   // output filename for display
    string       name;
    image<vec4f> img;
    bool         is_hdr = false;

    // diplay image
    image<vec4f> display;
    gltexture    gl_txt = {};

    // image stats
    image_stats stats;

    // tonemapping values
    float exposure = 0;
    bool  filmic   = false;
    bool  srgb     = true;

    // computation futures
    atomic<bool> load_done, display_done, stats_done, texture_done;
    thread       load_thread, display_thread, stats_thread, save_thread;
    atomic<bool> display_stop;
    string       error_msg = "";

    // viewing properties
    vec2f imcenter    = zero2f;
    float imscale     = 1;
    bool  zoom_to_fit = false;

    // cleanup and stop threads
    ~app_image() {
        display_stop = true;
        if (load_thread.joinable()) load_thread.join();
        if (display_thread.joinable()) display_thread.join();
        if (stats_thread.joinable()) stats_thread.join();
        if (save_thread.joinable()) save_thread.join();
    }
};

struct app_state {
    // images
    vector<app_image*> imgs;
    int                img_id = 0;

    ~app_state() {
        for (auto img : imgs) delete img;
    }
};

// compute min/max
void update_stats_async(app_image* img) {
    img->stats_done       = false;
    img->stats.pxl_bounds = invalid_bbox4f;
    img->stats.lum_bounds = invalid_bbox1f;
    for (auto p : img->img.pixels) {
        img->stats.pxl_bounds += p;
        img->stats.lum_bounds += luminance(xyz(p));
    }
    img->stats_done = true;
}

void update_display_async(app_image* img) {
    auto start        = get_time();
    img->display_done = false;
    img->texture_done = false;
    if (img->img.width * img->img.height > 1024 * 1024) {
        auto nthreads = thread::hardware_concurrency();
        auto threads  = vector<thread>();
        for (auto tid = 0; tid < nthreads; tid++) {
            threads.push_back(thread([img, tid, nthreads]() {
                for (auto j = tid; j < img->img.height; j += nthreads) {
                    if (img->display_stop) break;
                    for (auto i = 0; i < img->img.width; i++) {
                        pixel_at(img->display, i, j) = tonemap_filmic(
                            pixel_at(img->img, i, j), img->exposure,
                            img->filmic, img->srgb);
                    }
                }
            }));
        }
        for (auto& t : threads) t.join();
    } else {
        for (auto j = 0; j < img->img.height; j++) {
            if (img->display_stop) break;
            for (auto i = 0; i < img->img.width; i++) {
                pixel_at(img->display, i, j) = tonemap_filmic(
                    pixel_at(img->img, i, j), img->exposure, img->filmic,
                    img->srgb);
            }
        }
    }
    img->display_done = true;
    if (!img->display_stop) {
        printf("display: %s\n", format_duration(get_time() - start).c_str());
        fflush(stdout);
    }
}

// load image
void load_image_async(app_image* img) {
    auto start        = get_time();
    img->load_done    = false;
    img->stats_done   = false;
    img->display_done = false;
    img->texture_done = false;
    img->error_msg    = "";
    img->img          = load_image4f(img->filename);
    if (img->img.pixels.empty()) img->error_msg = "cannot load image";
    img->load_done = true;
    img->display   = img->img;
    printf("load: %s\n", format_duration(get_time() - start).c_str());
    fflush(stdout);
    img->display_thread = thread(update_display_async, img);
    img->stats_thread   = thread(update_stats_async, img);
}

// save an image
void save_image_async(app_image* img) {
    if (is_hdr_filename(img->outname)) {
        if (!save_image4b(img->outname, float_to_byte(img->display))) {
            img->error_msg = "error saving image";
        }
    } else {
        if (!save_image4f(img->outname, srgb_to_linear(img->display))) {
            img->error_msg = "error saving image";
        }
    }
}

// add a new image
void add_new_image(app_state* app, const string& filename, const string& outname,
    float exposure = 0, bool filmic = false, bool srgb = true) {
    auto img      = new app_image();
    img->filename = filename;
    img->outname = (outname == "") ? replace_extension(filename, ".display.png") :
                                     outname;
    img->name        = get_filename(filename);
    img->is_hdr      = is_hdr_filename(filename);
    img->exposure    = exposure;
    img->filmic      = filmic;
    img->srgb        = srgb;
    img->load_thread = thread(load_image_async, img);
    app->imgs.push_back(img);
    app->img_id = (int)app->imgs.size() - 1;
}

void draw_glwidgets(glwindow* win) {
    auto app    = (app_state*)get_user_pointer(win);
    auto edited = false;
    begin_glwidgets_frame(win);
    if (begin_glwidgets_window(win, "yimview")) {
        auto img = app->imgs.at(app->img_id);
        if (begin_header_glwidget(win, "image")) {
            draw_combobox_glwidget(win, "image", app->img_id, app->imgs);
            draw_label_glwidgets(win, "filename", "%s", img->filename.c_str());
            draw_textinput_glwidget(win, "outname", img->outname);
            if (draw_button_glwidget(win, "save display")) {
                if (img->display_done) {
                    img->save_thread = thread(save_image_async, img);
                }
            }
            auto status = string();
            if (img->error_msg != "")
                status = "error: " + img->error_msg;
            else if (!img->load_done)
                status = "loading...";
            else if (!img->display_done)
                status = "displaying...";
            else
                status = "done";
            draw_label_glwidgets(win, "status", status.c_str());
            draw_label_glwidgets(
                win, "size", "%d x %d ", img->img.width, img->img.height);
            draw_slider_glwidget(win, "zoom", img->imscale, 0.1, 10);
            draw_checkbox_glwidget(win, "zoom to fit", img->zoom_to_fit);
            end_header_glwidget(win);
        }
        if (begin_header_glwidget(win, "adjust")) {
            edited += draw_slider_glwidget(win, "exposure", img->exposure, -5, 5);
            edited += draw_checkbox_glwidget(win, "filmic", img->filmic);
            edited += draw_checkbox_glwidget(win, "srgb", img->srgb);
            end_header_glwidget(win);
        }
        if (begin_header_glwidget(win, "inspect")) {
            auto mouse_pos = get_glmouse_pos(win);
            auto ij = get_image_coords(mouse_pos, img->imcenter, img->imscale,
                {img->img.width, img->img.height});
            draw_dragger_glwidget(win, "mouse", ij);
            auto pixel = zero4f;
            if (ij.x >= 0 && ij.x < img->img.width && ij.y >= 0 &&
                ij.y < img->img.height) {
                pixel = pixel_at(img->img, ij.x, ij.y);
            }
            draw_coloredit_glwidget(win, "pixel", pixel);
            auto stats = (img->stats_done) ? img->stats : image_stats{};
            draw_dragger_glwidget(win, "pxl min", stats.pxl_bounds.min);
            draw_dragger_glwidget(win, "pxl max", stats.pxl_bounds.max);
            draw_dragger_glwidget(win, "lum min", stats.lum_bounds.min);
            draw_dragger_glwidget(win, "lum max", stats.lum_bounds.max);
            end_header_glwidget(win);
        }
    }
    end_glwidgets_frame(win);
    if (edited) {
        auto img = app->imgs.at(app->img_id);
        if (img->display_thread.joinable()) {
            img->display_stop = true;
            img->display_thread.join();
        }
        if (img->save_thread.joinable()) img->save_thread.join();
        img->display_stop   = false;
        img->display_thread = thread(update_display_async, img);
    }
}

void draw(glwindow* win) {
    auto app      = (app_state*)get_user_pointer(win);
    auto img      = app->imgs.at(app->img_id);
    auto win_size = get_glwindow_size(win);
    auto fb_size  = get_glframebuffer_size(win);
    set_glviewport(fb_size);
    clear_glframebuffer(vec4f{0.8f, 0.8f, 0.8f, 1.0f});
    if (img->gl_txt) {
        center_image4f(img->imcenter, img->imscale,
            {img->display.width, img->display.height}, win_size,
            img->zoom_to_fit);
        draw_glimage(img->gl_txt, {img->display.width, img->display.height},
            win_size, img->imcenter, img->imscale);
    }
    draw_glwidgets(win);
    swap_glbuffers(win);
}

void update(app_state* app) {
    for (auto img : app->imgs) {
        if (!img->display_done || img->texture_done) continue;
        if (!img->gl_txt) {
            img->gl_txt = make_gltexture(img->display, false, false, true);
        } else {
            update_gltexture(img->gl_txt, img->display, false, false, true);
        }
        img->texture_done = true;
    }
}

void drop_callback(glwindow* win, const vector<string>& paths) {
    auto app = (app_state*)get_user_pointer(win);
    for (auto path : paths) add_new_image(app, path, "");
}

void run_ui(app_state* app) {
    // window
    auto img   = app->imgs.at(app->img_id);
    auto width = 720 + 320, height = 720;
    auto win = make_glwindow(720 + 320, 720, "yimview", app, draw);
    set_drop_callback(win, drop_callback);

    // init widgets
    init_glwidgets(win);

    // center image
    center_image4f(img->imcenter, img->imscale, {img->img.width, img->img.height},
        {width, height}, img->img.width > width || img->img.height > height);

    // window values
    auto mouse_pos = zero2f, last_pos = zero2f;
    while (!should_glwindow_close(win)) {
        last_pos            = mouse_pos;
        mouse_pos           = get_glmouse_pos(win);
        auto mouse_left     = get_glmouse_left(win);
        auto mouse_right    = get_glmouse_right(win);
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
    auto parser   = make_cmdline_parser(argc, argv, "view images", "yimview");
    auto exposure = parse_arg(parser, "--exposure,-e", 0.0f, "display exposure");
    auto filmic   = parse_arg(parser, "--filmic", false, "display filmic");
    auto srgb     = parse_arg(parser, "--no-srgb", true, "display as sRGB");
    // auto quiet = parse_flag(
    //     parser, "--quiet,-q", false, "Print only errors messages");
    auto outfilename = parse_arg(parser, "--out,-o", ""s, "image out filename");
    auto filenames   = parse_args(
        parser, "images", vector<string>{}, "image filenames", true);
    check_cmdline(parser);

    // loading images
    for (auto filename : filenames)
        add_new_image(app, filename, outfilename, exposure, filmic, srgb);
    app->img_id = 0;

    // run ui
    run_ui(app);

    // cleanup
    delete app;

    // done
    return 0;
}
