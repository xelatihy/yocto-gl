//
// LICENSE:
//
// Copyright (c) 2016 -- 2019 Fabio Pellacini
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

#include "../yocto/yocto_image.h"
#include "../yocto/yocto_utils.h"
#include "yocto_opengl.h"
using namespace yocto;

#include "ext/CLI11.hpp"

struct image_stats {
    bbox4f        bounds    = {zero4f, zero4f};
    vec4f         average   = zero4f;
    vector<vec3f> histogram = {};
};

struct app_image {
    // original data
    string name     = "";
    string filename = "";
    string outname  = "";

    // image data
    image<vec4f> img = {};

    // diplay image
    image<vec4f>   display = {};
    opengl_texture gl_txt  = {};

    // image stats
    image_stats image_stats, display_stats;

    // tonemapping values
    tonemap_image_options    tonemap_options    = {};
    colorgrade_image_options colorgrade_options = {};

    // computation futures
    atomic<bool>                   load_done, display_done;
    thread                         load_thread, display_thread, save_thread;
    atomic<bool>                   display_stop;
    concurrent_queue<image_region> display_queue;
    string                         error_msg = "";

    // viewing properties
    vec2f image_center = zero2f;
    float image_scale  = 1;
    bool  zoom_to_fit  = false;

    // cleanup and stop threads
    ~app_image() {
        display_stop = true;
        if (load_thread.joinable()) load_thread.join();
        if (display_thread.joinable()) display_thread.join();
        if (save_thread.joinable()) save_thread.join();
    }
};

struct app_state {
    deque<app_image> images;
    int selected = -1;
};

// compute min/max
void compute_image_stats(
    image_stats& stats, const image<vec4f>& img, bool linear_hdr) {
    auto timer     = log_timed("computing stats");
    auto max_histo = linear_hdr ? 8 : 1;
    stats.bounds   = invalid_bbox4f;
    stats.average  = zero4f;
    stats.histogram.assign(256, zero3f);
    for (auto& p : img) {
        stats.bounds += p;
        stats.average += p;
        stats.histogram[(int)(clamp(p.x / max_histo, 0.f, 1.f) * 255)].x += 1;
        stats.histogram[(int)(clamp(p.y / max_histo, 0.f, 1.f) * 255)].y += 1;
        stats.histogram[(int)(clamp(p.z / max_histo, 0.f, 1.f) * 255)].z += 1;
    }
    auto num_pixels = (size_t)img.size().x * (size_t)img.size().y;
    for (auto& v : stats.histogram) v /= num_pixels;
    stats.average /= num_pixels;
}

void update_display_async(app_image& img) {
    auto timer       = log_timed("rendering {}", get_filename(img.filename));
    img.display_done = false;
    auto regions     = vector<image_region>{};
    make_image_regions(regions, img.img.size());
    parallel_foreach(
        regions,
        [&img,
            colorgrade = img.colorgrade_options != colorgrade_image_options{},
            tonemap_options    = img.tonemap_options,
            colorgrade_options = img.colorgrade_options](
            const image_region& region) {
            tonemap_image_region(img.display, region, img.img, tonemap_options);
            if (colorgrade) {
                colorgrade_image_region(
                    img.display, region, img.display, colorgrade_options);
            }
            img.display_queue.push(region);
        },
        &img.display_stop);
    compute_image_stats(img.display_stats, img.display, false);
    img.display_done = true;
}

// load image
void load_image_async(app_image& img) {
    auto timer       = log_timed("loading {}", get_filename(img.filename));
    img.name         = get_basename(img.filename) + " [loading]";
    img.load_done    = false;
    img.display_done = false;
    img.error_msg    = "";
    img.img          = {};
    try {
        load_image(img.filename, img.img);
        img.name = get_filename(img.filename) +
                   format(" [{}x{}]", img.img.size().x, img.img.size().y);
    } catch (const std::exception& e) {
        img.error_msg = e.what();
        img.name      = get_filename(img.filename) + " [error]";
        return;
    }
    compute_image_stats(
        img.image_stats, img.img, is_hdr_filename(img.filename));
    img.load_done      = true;
    img.display        = img.img;
    img.display_thread = thread([&img]() { update_display_async(img); });
    img.zoom_to_fit    = true;
}

// save an image
void save_image_async(app_image& img) {
    try {
        if (!is_hdr_filename(img.outname)) {
            auto ldr = image<vec4b>{};
            float_to_byte(ldr, img.display);
            save_image(img.outname, ldr);
        } else {
            auto aux = image<vec4f>{};
            srgb_to_linear(aux, img.display);
            save_image(img.outname, aux);
        }
    } catch (const std::exception& e) {
        img.error_msg = e.what();
    }
}

// add a new image
void add_new_image(app_state& app, const string& filename,
    const string& outname, const tonemap_image_options& tonemap_options = {}) {
    app.images.emplace_back();
    auto& img    = app.images.back();
    img.filename = filename;
    img.outname  = (outname == "") ? get_noextension(filename) + ".display.png"
                                  : outname;
    img.name            = get_filename(filename);
    img.tonemap_options = tonemap_options;
    if (!is_hdr_filename(filename)) {
        img.tonemap_options.filmic = false;
    }
    img.load_done    = false;
    img.display_done = false;
    app.selected     = (int)app.images.size() - 1;
}

void draw_opengl_widgets(const opengl_window& win) {
    auto& app    = *(app_state*)get_opengl_user_pointer(win);
    auto  edited = false;
    if (!begin_opengl_widgets_window(win, "yimview")) return;
    if (draw_button_opengl_widget(win, "load")) {
    }
    continue_opengl_widget_line(win);
    if (draw_button_opengl_widget(win, "save") && app.selected >= 0) {
        auto& img = app.images.at(app.selected);
        if (img.display_done) {
            img.save_thread = thread([&img]() { save_image_async(img); });
        }
    }
    continue_opengl_widget_line(win);
    if (draw_button_opengl_widget(win, "close") && app.selected >= 0) {
        // app.images.erase(app.images.begin()+app.selected);
        // app.selected = app.images.empty() ? -1 : 0;
    }
    continue_opengl_widget_line(win);
    if (draw_button_opengl_widget(win, "quit")) {
        set_close_opengl_window(win, true);
    }
    draw_combobox_opengl_widget(
        win, "image", app.selected, (int)app.images.size(),
        [&app](int idx) { return app.images[idx].name.c_str(); }, false);
    if (begin_header_opengl_widget(win, "tonemap")) {
        auto& img     = app.images.at(app.selected);
        auto  options = img.tonemap_options;
        draw_slider_opengl_widget(win, "exposure", options.exposure, -5, 5);
        draw_coloredit_opengl_widget(win, "tint", options.tint);
        draw_slider_opengl_widget(win, "contrast", options.contrast, 0, 1);
        draw_slider_opengl_widget(
            win, "logcontrast", options.logcontrast, 0, 1);
        draw_slider_opengl_widget(win, "saturation", options.saturation, 0, 1);
        draw_checkbox_opengl_widget(win, "filmic", options.filmic);
        continue_opengl_widget_line(win);
        draw_checkbox_opengl_widget(win, "srgb", options.srgb);
        continue_opengl_widget_line(win);
        if (draw_button_opengl_widget(win, "auto wb")) {
            edited       = true;
            auto wb      = 1 / xyz(img.image_stats.average);
            options.tint = wb / max(wb);
        }
        if (options != img.tonemap_options) {
            edited              = true;
            img.tonemap_options = options;
        }
        end_header_opengl_widget(win);
    }
    if (begin_header_opengl_widget(win, "colorgrade")) {
        auto& img     = app.images.at(app.selected);
        auto  options = img.colorgrade_options;
        draw_slider_opengl_widget(win, "contrast", options.contrast, 0, 1);
        draw_slider_opengl_widget(win, "ldr shadows", options.shadows, 0, 1);
        draw_slider_opengl_widget(win, "ldr midtones", options.midtones, 0, 1);
        draw_slider_opengl_widget(win, "highlights", options.highlights, 0, 1);
        draw_coloredit_opengl_widget(
            win, "shadows color", options.shadows_color);
        draw_coloredit_opengl_widget(
            win, "midtones color", options.midtones_color);
        draw_coloredit_opengl_widget(
            win, "highlights color", options.highlights_color);
        if (options != img.colorgrade_options) {
            edited                 = true;
            img.colorgrade_options = options;
        }
        end_header_opengl_widget(win);
    }
    if (begin_header_opengl_widget(win, "inspect")) {
        auto& img = app.images.at(app.selected);
        draw_label_opengl_widget(win, "filename", "%s", img.filename.c_str());
        draw_textinput_opengl_widget(win, "outname", img.outname);
        draw_slider_opengl_widget(win, "zoom", img.image_scale, 0.1, 10);
        draw_checkbox_opengl_widget(win, "zoom to fit", img.zoom_to_fit);
        auto mouse_pos = get_opengl_mouse_pos(win);
        auto ij        = get_image_coords(
            mouse_pos, img.image_center, img.image_scale, img.img.size());
        draw_dragger_opengl_widget(win, "mouse", ij);
        auto img_pixel = zero4f, display_pixel = zero4f;
        if (ij.x >= 0 && ij.x < img.img.size().x && ij.y >= 0 &&
            ij.y < img.img.size().y) {
            img_pixel     = img.img[{ij.x, ij.y}];
            display_pixel = img.display[{ij.x, ij.y}];
        }
        draw_coloredit_opengl_widget(win, "image", img_pixel);
        draw_dragger_opengl_widget(win, "display", display_pixel);
        auto img_stats = (img.load_done) ? img.image_stats : image_stats{};
        draw_dragger_opengl_widget(win, "image min", img_stats.bounds.min);
        draw_dragger_opengl_widget(win, "image max", img_stats.bounds.max);
        draw_dragger_opengl_widget(win, "image avg", img_stats.average);
        draw_histogram_opengl_widget(win, "image histo", img_stats.histogram);
        auto display_stats = (img.load_done) ? img.display_stats
                                             : image_stats{};
        draw_dragger_opengl_widget(
            win, "display min", display_stats.bounds.min);
        draw_dragger_opengl_widget(
            win, "display max", display_stats.bounds.max);
        draw_dragger_opengl_widget(win, "display avg", display_stats.average);
        draw_histogram_opengl_widget(
            win, "display histo", display_stats.histogram);
        end_header_opengl_widget(win);
    }
    if (begin_header_opengl_widget(win, "log")) {
        draw_log_opengl_widget(win);
    }
    if (edited) {
        auto& img = app.images.at(app.selected);
        if (img.display_thread.joinable()) {
            img.display_stop = true;
            img.display_thread.join();
        }
        if (img.save_thread.joinable()) img.save_thread.join();
        img.display_stop   = false;
        img.display_thread = thread([&img]() { update_display_async(img); });
    }
}

void draw(const opengl_window& win) {
    auto& app      = *(app_state*)get_opengl_user_pointer(win);
    auto  win_size = get_opengl_window_size(win);
    auto  fb_size  = get_opengl_framebuffer_size(win);
    set_opengl_viewport(fb_size);
    clear_opengl_lframebuffer(vec4f{0.15f, 0.15f, 0.15f, 1.0f});
    for(auto& img : app.images) {
        if (!img.load_done) continue;
        if (!img.gl_txt) {
            init_opengl_texture(
                img.gl_txt, img.img.size(), false, false, false, false);
        } else {
            auto region = image_region{};
            while (img.display_queue.try_pop(region)) {
                update_opengl_texture_region(
                    img.gl_txt, img.display, region, false);
            }
        }
    }
    if(!app.images.empty() && app.selected >= 0) {
        auto& img      = app.images.at(app.selected);
        if (img.load_done && img.gl_txt) {
            update_image_view(img.image_center, img.image_scale, img.display.size(),
                win_size, img.zoom_to_fit);
            draw_opengl_image_background(img.gl_txt, win_size.x, win_size.y,
                img.image_center, img.image_scale);
            set_opengl_blending(true);
            draw_opengl_image(img.gl_txt, win_size.x, win_size.y, img.image_center,
                img.image_scale);
            set_opengl_blending(false);
        }
    }
    begin_opengl_widgets_frame(win);
    draw_opengl_widgets(win);
    end_opengl_widgets_frame(win);
    swap_opengl_buffers(win);
}

void update(app_state& app) {}

void drop_callback(const opengl_window& win, const vector<string>& paths) {
    auto& app = *(app_state*)get_opengl_user_pointer(win);
    for (auto path : paths) {
        add_new_image(app, path, "");
        app.images.back().load_thread = thread(
            [&img = app.images.back()]() { load_image_async(img); });
    }
}

void run_ui(app_state& app) {
    // window
    auto win = opengl_window();
    init_opengl_window(win, {1280, 720}, "yimview", &app, draw);
    set_drop_opengl_callback(win, drop_callback);

    // init widgets
    init_opengl_widgets(win);

    // setup logging
    set_log_callback(
        [&win](const string& msg) { add_log_opengl_widget(win, msg.c_str()); });

    // load images
    for (auto& img : app.images) {
        img.load_thread = thread([&img]() { load_image_async(img); });
    }

    // window values
    auto mouse_pos = zero2f, last_pos = zero2f;
    while (!should_opengl_window_close(win)) {
        last_pos            = mouse_pos;
        mouse_pos           = get_opengl_mouse_pos(win);
        auto mouse_left     = get_opengl_mouse_left(win);
        auto mouse_right    = get_opengl_mouse_right(win);
        auto widgets_active = get_opengl_widgets_active(win);

        // handle mouse
        if (mouse_left && !widgets_active) {
            auto& img = app.images.at(app.selected);
            img.image_center += mouse_pos - last_pos;
        }
        if (mouse_right && !widgets_active) {
            auto& img = app.images.at(app.selected);
            img.image_scale *= powf(2, (mouse_pos.x - last_pos.x) * 0.001f);
        }

        // update
        update(app);

        // draw
        draw(win);

        // event hadling
        // could also wait if (mouse_left || mouse_right || widgets_active)
        process_opengl_events(win);
    }

    // cleanup
    delete_opengl_window(win);
}

int main(int argc, char* argv[]) {
    // prepare application
    auto app             = app_state();
    auto tonemap_options = tonemap_image_options{};
    auto outfilename     = ""s;
    auto filenames       = vector<string>{};

    // command line options
    auto parser = CLI::App{"view images"};
    parser.add_option(
        "--exposure,-e", tonemap_options.exposure, "display exposure");
    parser.add_flag(
        "--filmic,!--no-filmic", tonemap_options.filmic, "display filmic");
    parser.add_flag(
        "--srgb,!--no-srgb", tonemap_options.srgb, "display as sRGB");
    // auto quiet = parse_flag(
    //     parser, "--quiet,-q", false, "Print only errors messages");
    parser.add_option("--out,-o", outfilename, "image out filename");
    parser.add_option("images", filenames, "image filenames")->required(true);
    try {
        parser.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        return parser.exit(e);
    }

    // loading images
    for (auto filename : filenames)
        add_new_image(app, filename, outfilename, tonemap_options);
    app.selected = 0;

    // run ui
    run_ui(app);

    // done
    return 0;
}
