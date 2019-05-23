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
    vec4f         min       = zero4f;
    vec4f         max       = zero4f;
    vec4f         average   = zero4f;
    vector<vec3f> histogram = {};
};

enum struct app_task_type { none, load, save, display, close };

struct app_task {
    app_task_type                  type;
    future<void>                   result;
    atomic<bool>                   stop;
    concurrent_queue<image_region> queue;

    app_task(app_task_type type) : type{type}, result{}, stop{false}, queue{} {}
    ~app_task() {
        stop = true;
        if (result.valid()) {
            try {
                result.get();
            } catch (...) {
            }
        }
    }
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
    tonemap_image_params    tonemap_params    = {};
    colorgrade_image_params colorgrade_params = {};

    // computation futures
    bool            load_done = false, display_done = false;
    deque<app_task> task_queue;

    // viewing properties
    vec2f image_center = zero2f;
    float image_scale  = 1;
    bool  zoom_to_fit  = false;
};

struct app_state {
    // data
    deque<app_image> images;
    int              selected = -1;
    deque<string>    errors;

    // default options
    tonemap_image_params    tonemap_params    = {};
    colorgrade_image_params colorgrade_params = {};
};

// compute min/max
void compute_image_stats(
    image_stats& stats, const image<vec4f>& img, bool linear_hdr) {
    auto max_histo = linear_hdr ? 8 : 1;
    stats.min      = vec4f{float_max};
    stats.max      = vec4f{float_min};
    stats.average  = zero4f;
    stats.histogram.assign(256, zero3f);
    for (auto& p : img) {
        stats.min = min(stats.min, p);
        stats.max = max(stats.max, p);
        stats.average += p;
        stats.histogram[(int)(clamp(p.x / max_histo, 0.f, 1.f) * 255)].x += 1;
        stats.histogram[(int)(clamp(p.y / max_histo, 0.f, 1.f) * 255)].y += 1;
        stats.histogram[(int)(clamp(p.z / max_histo, 0.f, 1.f) * 255)].z += 1;
    }
    auto num_pixels = (size_t)img.size().x * (size_t)img.size().y;
    for (auto& v : stats.histogram) v /= num_pixels;
    stats.average /= num_pixels;
}

void update_app_display(const string& filename, const image<vec4f>& img,
    image<vec4f>& display, image_stats& stats,
    const tonemap_image_params&    tonemap_params,
    const colorgrade_image_params& colorgrade_params, atomic<bool>& stop,
    concurrent_queue<image_region>& queue) {
    auto regions = vector<image_region>{};
    make_regions(regions, img.size(), 128);
    parallel_foreach(
        regions,
        [&img, &display, &queue, tonemap_params, colorgrade_params,
            do_colorgrade = colorgrade_params != colorgrade_image_params{}](
            const image_region& region) {
            tonemap(display, img, region, tonemap_params);
            if (do_colorgrade) {
                colorgrade(display, display, region, colorgrade_params);
            }
            queue.push(region);
        },
        &stop);
    compute_image_stats(stats, display, false);
}

// add a new image
void add_new_image(app_state& app, const string& filename) {
    auto& img             = app.images.emplace_back();
    img.filename          = filename;
    img.outname           = get_noextension(filename) + ".display.png";
    img.name              = get_filename(filename);
    img.tonemap_params    = app.tonemap_params;
    img.colorgrade_params = app.colorgrade_params;
    img.load_done         = false;
    img.display_done      = false;
    img.task_queue.emplace_back(app_task_type::load);
    app.selected = (int)app.images.size() - 1;
}

void draw_opengl_widgets(const opengl_window& win) {
    static string load_path = "", save_path = "", error_message = "";
    auto&         app = *(app_state*)get_opengl_user_pointer(win);
    if (!begin_opengl_widgets_window(win, "yimview")) return;
    if (!app.errors.empty() && error_message.empty()) {
        error_message = app.errors.front();
        app.errors.pop_front();
        open_modal_opengl_widget(win, "error");
    }
    if (!draw_modal_message_opengl_window(win, "error", error_message)) {
        error_message = "";
    }
    if (draw_modal_fileialog_opengl_widgets(win, "load image", load_path, false,
            "./", "", "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
        add_new_image(app, load_path);
    }
    if (draw_modal_fileialog_opengl_widgets(win, "save image", save_path, true,
            get_dirname(save_path), get_filename(save_path),
            "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
        app.images[app.selected].outname = save_path;
        app.images[app.selected].task_queue.emplace_back(app_task_type::save);
        save_path = "";
    }
    if (draw_button_opengl_widget(win, "load")) {
        open_modal_opengl_widget(win, "load image");
    }
    continue_opengl_widget_line(win);
    if (draw_button_opengl_widget(win, "save",
            app.selected >= 0 && app.images[app.selected].display_done)) {
        save_path = app.images[app.selected].outname;
        open_modal_opengl_widget(win, "save image");
    }
    continue_opengl_widget_line(win);
    if (draw_button_opengl_widget(win, "close", app.selected >= 0)) {
        auto& img = app.images.at(app.selected);
        img.task_queue.emplace_back(app_task_type::close);
    }
    continue_opengl_widget_line(win);
    if (draw_button_opengl_widget(win, "quit")) {
        set_close_opengl_window(win, true);
    }
    if (app.images.empty()) return;
    draw_combobox_opengl_widget(
        win, "image", app.selected, (int)app.images.size(),
        [&app](int idx) { return app.images[idx].name.c_str(); }, false);
    auto& img = app.images.at(app.selected);
    if (begin_header_opengl_widget(win, "tonemap")) {
        auto options = img.tonemap_params;
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
            auto wb      = 1 / xyz(img.image_stats.average);
            options.tint = wb / max(wb);
        }
        if (options != img.tonemap_params) {
            img.tonemap_params = options;
            if (img.load_done)
                img.task_queue.emplace_back(app_task_type::display);
        }
        end_header_opengl_widget(win);
    }
    if (begin_header_opengl_widget(win, "colorgrade")) {
        auto options = img.colorgrade_params;
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
        if (options != img.colorgrade_params) {
            img.colorgrade_params = options;
            if (img.load_done)
                img.task_queue.emplace_back(app_task_type::display);
        }
        end_header_opengl_widget(win);
    }
    if (begin_header_opengl_widget(win, "inspect")) {
        draw_label_opengl_widget(win, "image", get_filename(img.filename));
        draw_label_opengl_widget(win, "filename", img.filename);
        draw_label_opengl_widget(win, "outname", img.outname);
        draw_label_opengl_widget(
            win, "image", "%d x %d", img.img.size().x, img.img.size().y);
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
        draw_dragger_opengl_widget(win, "image min", img_stats.min);
        draw_dragger_opengl_widget(win, "image max", img_stats.max);
        draw_dragger_opengl_widget(win, "image avg", img_stats.average);
        draw_histogram_opengl_widget(win, "image histo", img_stats.histogram);
        auto display_stats = (img.load_done) ? img.display_stats
                                             : image_stats{};
        draw_dragger_opengl_widget(win, "display min", display_stats.min);
        draw_dragger_opengl_widget(win, "display max", display_stats.max);
        draw_dragger_opengl_widget(win, "display avg", display_stats.average);
        draw_histogram_opengl_widget(
            win, "display histo", display_stats.histogram);
        end_header_opengl_widget(win);
    }
    if (begin_header_opengl_widget(win, "log")) {
        draw_log_opengl_widget(win);
        end_header_opengl_widget(win);
    }
}

void draw(const opengl_window& win) {
    auto& app      = *(app_state*)get_opengl_user_pointer(win);
    auto  win_size = get_opengl_window_size(win);
    auto  fb_view  = get_opengl_framebuffer_viewport(win);
    set_opengl_viewport(fb_view);
    clear_opengl_framebuffer(vec4f{0.15f, 0.15f, 0.15f, 1.0f});
    if (!app.images.empty() && app.selected >= 0) {
        auto& img = app.images.at(app.selected);
        if (img.load_done && img.gl_txt) {
            update_image_view(img.image_center, img.image_scale,
                img.display.size(), win_size, img.zoom_to_fit);
            draw_opengl_image_background(img.gl_txt, win_size.x, win_size.y,
                img.image_center, img.image_scale);
            set_opengl_blending(true);
            draw_opengl_image(img.gl_txt, win_size.x, win_size.y,
                img.image_center, img.image_scale);
            set_opengl_blending(false);
        }
    }
    begin_opengl_widgets_frame(win);
    draw_opengl_widgets(win);
    end_opengl_widgets_frame(win);
    swap_opengl_buffers(win);
}

void update(app_state& app) {
    // close if needed
    while (!app.images.empty()) {
        auto pos = -1;
        for (auto idx = 0; idx < app.images.size(); idx++) {
            for (auto& task : app.images[idx].task_queue) {
                if (task.type == app_task_type::close) pos = idx;
            }
        }
        if (pos < 0) break;
        app.images.erase(app.images.begin() + pos);
        app.selected = app.images.empty() ? -1 : 0;
    }
    // consume partial results
    for (auto& img : app.images) {
        if (img.task_queue.empty()) continue;
        auto& task = img.task_queue.front();
        if (task.type != app_task_type::display || task.queue.empty()) continue;
        auto region = image_region{};
        while (img.task_queue.front().queue.try_pop(region)) {
            update_opengl_texture_region(
                img.gl_txt, img.display, region, false);
        }
    }
    // remove unneeded tasks
    for (auto& img : app.images) {
        while (img.task_queue.size() > 1) {
            auto& task = img.task_queue.at(0);
            auto& next = img.task_queue.at(1);
            if (task.type == app_task_type::display) {
                if (next.type != app_task_type::display) break;
                log_info("cancel rendering " + img.filename);
            } else {
                break;
            }
            task.stop = true;
            if (task.result.valid()) {
                try {
                    task.result.get();
                } catch (...) {
                }
            }
            img.task_queue.pop_front();
        }
    }
    // grab result of finished tasks
    for (auto& img : app.images) {
        if (img.task_queue.empty()) continue;
        auto& task = img.task_queue.front();
        if (!task.result.valid() || !task.queue.empty() ||
            task.result.wait_for(std::chrono::nanoseconds(10)) !=
                std::future_status::ready)
            continue;
        switch (task.type) {
            case app_task_type::none: break;
            case app_task_type::close: break;
            case app_task_type::load: {
                try {
                    task.result.get();
                    img.load_done = true;
                    img.name      = get_filename(img.filename) + " [" +
                               to_string(img.img.size()) + "]";
                    img.display = img.img;
                    log_info("done loading " + img.filename);
                    init_opengl_texture(
                        img.gl_txt, img.display, false, false, false);
                    img.task_queue.emplace_back(app_task_type::display);
                } catch (std::exception& e) {
                    log_error(e.what());
                    img.name = get_filename(img.filename) + " [error]";
                    app.errors.push_back("cannot load " + img.filename);
                }
            } break;
            case app_task_type::save: {
                try {
                    task.result.get();
                    log_info("done saving " + img.outname);
                } catch (std::exception& e) {
                    log_error(e.what());
                    app.errors.push_back("cannot save " + img.outname);
                }
            } break;
            case app_task_type::display: {
                try {
                    task.result.get();
                    img.display_done = true;
                    log_info("done rendering " + img.filename);
                } catch (std::exception& e) {
                    log_error(e.what());
                    app.errors.push_back("cannot render " + img.filename);
                }
            } break;
        }
        img.task_queue.pop_front();
    }
    // schedule tasks not running
    for (auto& img : app.images) {
        if (img.task_queue.empty()) continue;
        auto& task = img.task_queue.front();
        if (task.result.valid()) continue;
        task.stop = false;
        switch (task.type) {
            case app_task_type::none: break;
            case app_task_type::close: break;
            case app_task_type::load: {
                log_info("start loading " + img.filename);
                img.load_done = false;
                task.result   = async([&img]() {
                    img.img = {};
                    load_image(img.filename, img.img);
                    compute_image_stats(img.image_stats, img.img,
                        is_hdr_filename(img.filename));
                });
            } break;
            case app_task_type::save: {
                log_info("start saving " + img.outname);
                task.result = async([&img]() {
                    if (!is_hdr_filename(img.outname)) {
                        auto ldr = image<vec4b>{};
                        float_to_byte(ldr, img.display);
                        save_image(img.outname, ldr);
                    } else {
                        auto aux = image<vec4f>{};
                        srgb_to_rgb(aux, img.display);
                        save_image(img.outname, aux);
                    }
                });
            } break;
            case app_task_type::display: {
                log_info("start rendering " + img.filename);
                img.display_done = false;
                task.result      = async([&img, &task]() {
                    update_app_display(img.filename, img.img, img.display,
                        img.display_stats, img.tonemap_params,
                        img.colorgrade_params, task.stop, task.queue);
                });
            } break;
        }
    }
}

void drop_callback(const opengl_window& win, const vector<string>& paths) {
    auto& app = *(app_state*)get_opengl_user_pointer(win);
    for (auto path : paths) add_new_image(app, path);
}

void run_ui(app_state& app) {
    // window
    auto win = opengl_window();
    init_opengl_window(win, {1280 + 320, 720}, "yimview", &app, draw);
    set_drop_opengl_callback(win, drop_callback);

    // init widgets
    init_opengl_widgets(win);

    // setup logging
    set_log_callback(
        [&win](const string& str) { add_log_opengl_widget(win, str.c_str()); });

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
        process_opengl_events(win);
    }

    // cleanup
    delete_opengl_window(win);
}

int main(int argc, char* argv[]) {
    // prepare application
    auto app       = app_state();
    auto filenames = vector<string>{};

    // command line options
    auto parser = CLI::App{"view images"};
    // auto quiet = parse_flag(
    //     parser, "--quiet,-q", false, "Print only errors messages");
    parser.add_option("images", filenames, "image filenames")->required(true);
    try {
        parser.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        return parser.exit(e);
    }

    // loading images
    for (auto filename : filenames) add_new_image(app, filename);
    app.selected = 0;

    // run ui
    run_ui(app);

    // done
    return 0;
}
